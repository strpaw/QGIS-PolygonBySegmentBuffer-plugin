# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PolygonBySegmentBuffer
                                 A QGIS plugin
 Creates polygon which is defined as buffer around segment
                              -------------------
        begin                : 2018-08-08
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Paweł Strzelewicz
        email                : pawel.strzelewicz83@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt4.QtCore import *
from PyQt4.QtGui import QAction, QIcon, QFileDialog, QMessageBox, QWidget
from qgis.core import *
from qgis.gui import *
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from polygon_by_segmentbuffer_dialog import PolygonBySegmentBufferDialog
import os.path
import re
import math
import datetime

# Constants

# Parameters of WGS84 ellipsoid
WGS84_A = 6378137.0  # semi-major axis of the WGS84 ellipsoid in m
WGS84_B = 6356752.314245  # semi-minor axis of the WGS84 ellipsoid in m
WGS84_F = 1 / 298.257223563  # flattening of the WGS84 ellipsoid

# Special constants to use instead of False, to avoid ambigous where result of function might equal 0 and
# and result of function will be used in if statements etc.
VALID = 'VALID'
NOT_VALID = 'NOT_VALID'

# Units of measure
UOM_M = 'M'     # meters
UOM_KM = 'KM'    # kilometers
UOM_NM = 'NM'    # nautical miles
UOM_FEET = 'FEET'  # feet
UOM_SM = 'SM'    # statue miles

# Value types constant
V_AZM = 'AZM'
V_MAG_VAR = 'MAG_VAR'
V_LAT = 'LAT'
V_LON = 'LON'

# DMS, DM format separators, e.g. N32-44-55.21, N32 44 55.21, N32DEG 44 55.21
S_SPACE = ' '     # Blank, space separator
S_HYPHEN = '-'     # Hyphen separator
S_WORD_DEG = 'DEG'
S_WORD_MIN = 'MIN'
S_WORD_SEC = 'SEC'
S_ALL = [S_SPACE, S_HYPHEN, S_WORD_DEG, S_WORD_MIN, S_WORD_SEC]

# Hemisphere letters
H_ALL = ['N', 'S', 'E', 'W']
H_LAT = ['N', 'S']
H_LON = ['E', 'W']
H_MINUS = ['S', 'W']

def get_tmp_name():
    """ Creates temporary name based on current time """
    c_time = datetime.datetime.now()               # Current time
    tmp_name = str(c_time).replace('-', '')         # Remove hyphens
    tmp_name = tmp_name.replace(':','')             # Remove colons
    tmp_name = tmp_name.replace(' ', '_')
    tmp_name = tmp_name.replace('.', '_')
    tmp_name = tmp_name[2:] # Trim two first digits of year
    return tmp_name


def check_range_azimuth(azm):
    if azm < 0 or azm > 360:
        result = NOT_VALID
    else:
        result = azm
    return result


def check_azimuth_DD(azm):
    """ Azimuth in DD (decimal degrees) format validation.
    :param azm: string, azimuth to validate
    :return result: float, azimuth value in decima degress format if azm valu is valid azimuth,
                    constant NOT_VALID if azm value is not valid azimuth
    """
    try:
        a = float(azm)
        result = check_range_azimuth(a)
    except:
        result = NOT_VALID
    return result

def check_azimuth(azm):
    """ Azimuth in various format validation.
    :param azm: string, azimuth to validate
    :return result: float, azimuth value in decimal degress format if azm value is valid azimuth,
                    constant NOT_VALID if azm value is not valid azimuth
    """
    azm = str(azm)
    azm.strip()
    azm = azm.replace(",", ".")
    result = check_azimuth_DD(azm)
    return result



def check_range(value, value_type):
    """ Check if given value of given value_type is within the range for this value type.
    :param value: float, value to check range
    :param value_type : constant of value type (e.g V_AZM)
    """

    if value_type == V_LAT:
        if (value < -90) or (value > 90):
            result = NOT_VALID
        else:
            result = value
    elif value_type == V_LON:
        if value < -180 or value > 180:
            result = NOT_VALID
        else:
            result = value
    return result


def check_signed_dd(c, c_type):
    """ Checks if input parameter latitude or longitude in signed decimal degrees format
    :param c: string,
    :param c_type: type of coordinate (latitude or longitude)
    return dd: float - decimal degrees if is float value, NOT_VALID constant if is not valid float value
    """
    try:
        d = float(c)
        dd = check_range(d, c_type)
    except:
        dd = NOT_VALID
    return dd


def check_hletter_delimited_dms_dm(c, c_type):
    """ Checks if input parameter is DMS (degrees, minutes, seconds) or DM (degrees, minutes) with hemisphere letter prefix or suffix,
    :param c: string, coordinate to check
    :param c_type: string: type of coordinate (latitude or longitude)
    :return dd: float - decimal degrees if is valid dms, NOT_VALID constant if is not valid float value
    """
    dd = ''
    p_hem = c[0]
    s_hem = c[len(c) - 1]
    if p_hem not in H_ALL and s_hem not in H_ALL:
        dd = NOT_VALID
    elif p_hem in H_ALL and s_hem in H_ALL:
        dd = NOT_VALID
    elif p_hem in H_ALL and s_hem not in H_ALL:
        dms_n = c[1:] + p_hem
    else:
        dms_n = c

    if dd != NOT_VALID:
        h = dms_n[len(dms_n) - 1]
        dms_n = dms_n[0:len(dms_n) - 1]  # Trim hemisphere letter

        for sep in S_ALL:  # Replace any separator to space separator
            dms_n = dms_n.replace(sep, ' ')

        dms_n = dms_n.strip()  # Trim trailing spaces - when seconds sign ('') replaced by space at the end of the string
        dms_d = re.sub(r'\s+', ' ', dms_n)  # Replace multiple spaces to single space
        dms_t = dms_d.split(' ')  # Splits dms by spaces and return as tuple

        if len(dms_t) == 3:  # 3 elements in tuple - assumes it is DMS format (DD MM SS.sss)
            try:
                d = float(dms_t[0])
                m = float(dms_t[1])
                s = float(dms_t[2])
                if d < 0 or m < 0 or m >= 60 or s < 0 or s >= 60:
                    dd = NOT_VALID
                elif c_type == V_LAT and h not in H_LAT:
                    dd = NOT_VALID
                elif c_type == V_LON and h not in H_LON:
                    dd = NOT_VALID
                else:
                    dd = d + m / 60 + s / 3600
                    dd = check_range(dd, c_type)
                    if (h in H_MINUS) and (dd != NOT_VALID):
                        dd = -dd
            except:
                dd = NOT_VALID
        elif len(dms_t) == 2:  # 2 elements in tuple - assumes it is DM format (DD MM.mmmm)
            try:
                d = float(dms_t[0])
                m = float(dms_t[1])
                if (d < 0) or (m < 0) or (m >= 60):
                    dd = NOT_VALID
                elif c_type == V_LAT and h not in H_LAT:
                    dd = NOT_VALID
                elif c_type == V_LON and h not in H_LON:
                    dd = NOT_VALID
                else:
                    dd = d + m / 60
                    dd = check_range(dd, c_type)
                    if (h in H_MINUS) and (dd != NOT_VALID):
                        dd = -dd
            except:
                dd = NOT_VALID
        else:
            dd = NOT_VALID
    return dd


def parse_dms2dd(c, c_type):
    """ Checks if input parameter is float number, doesn't check latitude, longitude limiest (-90 +90, -180 +180)
    :param c: string
    :param c_type: type of coordinate (latitude or longitude)
    return dd: decimal degrees or NOT_VALID constant if input is not valid coordinate,
    """
    dms = str(c)  # Ensure that dms in string variable to perform some string built-in functions

    if dms == '':  # Empty string
        dd = NOT_VALID
    else:
        dms = dms.replace(',', '.')  # Replace comma decimal separator to period decimal separator
        dms = dms.strip(' ')  # Remove leading and trailing blanks
        dms = dms.upper()  # Ensure that all characters are in upper case, e.g N, E instead of n, e
        dd = check_signed_dd(dms, c_type)
        if dd == NOT_VALID:
            dd = check_hletter_delimited_dms_dm(dms, c_type)
    return dd


# Distance functions
# Pattern for distance regular expression
REGEX_DIST = re.compile(r'^\d+(\.\d+)?$')


def check_distance(d):
    """ Distance validation.
    :param d: string, distance to validate
    :return is_valid: constant VALID if distance is valid, constant NOT_VALID if distance is not valid (e.g distance is less than 0)
    """
    if REGEX_DIST.match(d):
        result = VALID
    else:
        result = NOT_VALID
    return result


def km2m(km):
    """ Converts kilometers to meters
    :param km: float, value in kilometers
    :return: value in meters
    """
    return km * 1000


def NM2m(NM):
    """ Converts nautical miles to meters
    :param NM: float, value in nautical miles
    :return: value in meters
    """
    return NM * 1852


def feet2m(feet):
    """ Converts feet to meters
    :param feet: float, value in feet
    :return: value in meters
    """
    return feet * 0.3048


def SM2m(sm):
    """ Converts statue miles to meters
    :param sm: float, value in statue miles
    :return: value in meters
    """
    return sm * 1609.344


def to_meters(d, d_unit):
    """ Converts distance given in feet, nautical miles, statue miles etc. to distance in meters
    :param d: float, distance
    :param d_unit: constant unit of measure, unit of measure
    :return dm: float, distance in meters
    """
    if d_unit == UOM_M:
        dm = d
    elif d_unit == UOM_KM:
        dm = d * 1000
    elif d_unit == UOM_FEET:
        dm = feet2m(d)
    elif d_unit == UOM_SM:
        dm = SM2m(d)
    elif d_unit == UOM_NM:
        dm = NM2m(d)
    return dm


def vincenty_direct_solution(begin_lat, begin_lon, begin_azimuth, distance, a, b, f):
    """ Computes the latitude and longitude of the second point based on latitude, longitude,
    of the first point and distance and azimuth from first point to second point.
    Uses the algorithm by Thaddeus Vincenty for direct geodetic problem.
    For more information refer to: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf

    :param begin_lat: float, latitude of the first point; decimal degrees
    :param begin_lon: float, longitude of the first point; decimal degrees
    :param begin_azimuth: float, azimuth from first point to second point; decimal degrees
    :param distance: float, distance from first point to second point; meters
    :param a: float, semi-major axis of ellipsoid; meters
    :param b: float, semi-minor axis of ellipsoid; meters
    :param f: float, flattening of ellipsoid
    :return lat2_dd, lon2_dd: float, float latitude and longitude of the second point, decimal degrees
    """
    # Convert latitude, longitude, azimuth of the begining point to radians
    lat1 = math.radians(begin_lat)
    lon1 = math.radians(begin_lon)
    alfa1 = math.radians(begin_azimuth)

    sinAlfa1 = math.sin(alfa1)
    cosAlfa1 = math.cos(alfa1)

    # U1 - reduced latitude
    tanU1 = (1 - f) * math.tan(lat1)
    cosU1 = 1 / math.sqrt(1 + tanU1 * tanU1)
    sinU1 = tanU1 * cosU1

    # sigma1 - angular distance on the sphere from the equator to begining point
    sigma1 = math.atan2(tanU1, math.cos(alfa1))

    # sinAlfa - azimuth of the geodesic at the equator
    sinAlfa = cosU1 * sinAlfa1
    cosSqAlfa = 1 - sinAlfa * sinAlfa
    uSq = cosSqAlfa * (a * a - b * b) / (b * b)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))

    sigma = distance / (b * A)
    sigmap = 1

    while (math.fabs(sigma - sigmap) > 1e-12):
        cos2sigmaM = math.cos(2 * sigma1 + sigma)
        sinSigma = math.sin(sigma)
        cosSigma = math.cos(sigma)
        dSigma = B * sinSigma * (cos2sigmaM + B / 4 * (
                    cosSigma * (-1 + 2 * cos2sigmaM * cos2sigmaM) - B / 6 * cos2sigmaM * (
                        -3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2sigmaM * cos2sigmaM)))
        sigmap = sigma
        sigma = distance / (b * A) + dSigma

    var_aux = sinU1 * sinSigma - cosU1 * cosSigma * cosAlfa1  # Auxiliary variable

    # Latitude of the end point in radians
    lat2 = math.atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlfa1,
                      (1 - f) * math.sqrt(sinAlfa * sinAlfa + var_aux * var_aux))

    lamb = math.atan2(sinSigma * sinAlfa1, cosU1 * cosSigma - sinU1 * sinSigma * cosAlfa1)
    C = f / 16 * cosSqAlfa * (4 + f * (4 - 3 * cosSqAlfa))
    L = lamb - (1 - C) * f * sinAlfa * (
                sigma + C * sinSigma * (cos2sigmaM + C * cosSigma * (-1 + 2 * cos2sigmaM * cos2sigmaM)))
    # Longitude of the second point in radians
    lon2 = (lon1 + L + 3 * math.pi) % (2 * math.pi) - math.pi

    # Convert to decimal degrees
    lat2_dd = math.degrees(lat2)
    lon2_dd = math.degrees(lon2)

    return lat2_dd, lon2_dd


def dist_azm_orth_offset2latlon(ref_lat, ref_lon, ref_azm, distance_m, offset_m, offset_side):
    """ Calculates latitude and longitude of the second point base don latitude, longitude of the firts point, azimuth, distance and orthogonal offset
    Example: distance 1500 m, azimuth 45 degress and offset 500 meter left
    :param ref_lat: float, reference point latitude
    :param ref_lon: float, reference poitn longitude
    :param ref_azm: float, azimuth from reference point to intermediate point
    :param distance_m: float, distance in meters
    :param offset_m: float, offset in meters
    :param offset_side: indicate offset side, 'LEFT' for left, 'RIGHT' for right
    :return lat2_dd, lon2_dd: float, second point latitude, longitude
    """
    # Calculate azimuth from intermediate point to second point
    if offset_side == 'LEFT':
        offset_azm = ref_azm - 90
    elif offset_side == 'RIGHT':
        offset_azm = ref_azm + 90
    # Normalize azm to [0,360] degrees
    if offset_azm < 0:
        offset_azm += 360
    elif offset_azm > 360:
        offset_azm -= 360

    # Calculate intermediate point latitude, longitude
    inter_lat_dd, inter_lon_dd = vincenty_direct_solution(ref_lat, ref_lon, ref_azm, distance_m, WGS84_A, WGS84_B,
                                                          WGS84_F)

    # Calculate second point latitude, longitude, as reference point use intermediate point
    lat2_dd, lon2_dd = vincenty_direct_solution(inter_lat_dd, inter_lon_dd, offset_azm, offset_m, WGS84_A, WGS84_B,
                                                WGS84_F)

    return lat2_dd, lon2_dd

w = QWidget()


class PolygonBySegmentBuffer:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgisInterface
        """

        self.polygon_name = ""
        self.ref_name = ""
        self.ref_lat = None
        self.ref_lon = None
        self.ref_magvar = None
        self.brng = None
        self.dist_a = None
        self.dist_b = None
        self.buffer = None
        self.mlyr_name = ""

        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'PolygonBySegmentBuffer_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)


        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&PolygonBySegmentBuffer')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'PolygonBySegmentBuffer')
        self.toolbar.setObjectName(u'PolygonBySegmentBuffer')

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('PolygonBySegmentBuffer', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        # Create the dialog (after translation) and keep reference
        self.dlg = PolygonBySegmentBufferDialog()

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/PolygonBySegmentBuffer/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'PolygonBySegmentBuffer'),
            callback=self.run,
            parent=self.iface.mainWindow())

        self.dlg.pbCreateRectangleBuffer.clicked.connect(self.create_rectangle_buffer)

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&PolygonBySegmentBuffer'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def check_input(self):
        check_result = True
        err_msg = ""

        if self.dlg.leRefLat.text() == "":
            check_result = False
            err_msg += "Enter latitude of reference point!\n"

        if self.dlg.leRefLat.text() != "":
            self.ref_lat = parse_dms2dd(self.dlg.leRefLat.text(), V_LAT)
            if self.ref_lat == NOT_VALID:
                check_result = False
                err_msg += 'Latitude latitude of reference point wrong format!\n'

        if self.dlg.leRefLon.text() == "":
            check_result = False
            err_msg += "Enter longitude of reference point!\n"

        if self.dlg.leRefLon.text() != "":
            self.ref_lon = parse_dms2dd(self.dlg.leRefLon.text(), V_LON)
            if self.ref_lon == NOT_VALID:
                check_result = False
                err_msg += 'Longitude longitude of reference point wrong format!\n'

        if self.dlg.leBrng.text() == '':
            err_msg += 'Enter azimuth/bearing!\n'
            check_result = False
        else:
            self.brng = check_azimuth(self.dlg.leBrng.text())
            if self.brng == NOT_VALID:
                err_msg += 'Azimuth/Bearing wrong format!\n'
                check_result = False

        if self.dlg.leDistA.text() == '':
            err_msg += 'Enter distance A!\n'
            check_result = False
        else:
            if check_distance(self.dlg.leDistA.text()) == NOT_VALID:
                err_msg += 'Distance A wrong format!\n'
                check_result = False
            else:
                self.dist_a = float(self.dlg.leDistA.text())

        if self.dlg.leDistB.text() == '':
            err_msg += 'Enter distance B!\n'
            check_result = False
        else:
            if check_distance(self.dlg.leDistB.text()) == NOT_VALID:
                err_msg += 'Distance B wrong format!\n'
                check_result = False
            else:
                self.dist_b = float(self.dlg.leDistB.text())

        if self.dlg.leBuffer.text() == '':
            err_msg += 'Enter buffer distance!\n'
            check_result = False
        else:
            if check_distance(self.dlg.leBuffer.text()) == NOT_VALID:
                err_msg += 'Distance B wrong format!\n'
                check_result = False
            else:
                self.buffer = float(self.dlg.leBuffer.text())


        if not check_result:
            QMessageBox.critical(w, "Message", err_msg)

        return check_result

    def get_buffer_uom(self):
        """ Get buffer unit of measure """
        if self.dlg.cbBufferUOM.currentIndex() == 0:  # m
            buffer_uom = UOM_M
        elif self.dlg.cbBufferUOM.currentIndex() == 1:  # KM
            buffer_uom = UOM_KM
        elif self.dlg.cbBufferUOM.currentIndex() == 2:  # NM
            buffer_uom = UOM_NM
        elif self.dlg.cbBufferUOM.currentIndex() == 3:  # feet
            buffer_uom = UOM_FEET
        elif self.dlg.cbBufferUOM.currentIndex() == 4:  # SM
            buffer_uom = UOM_SM
        return buffer_uom

    def get_dist_uom(self):
        """ Get distance unit of measure """
        if self.dlg.cbDistUOM.currentIndex() == 0:  # m
            dist_uom = UOM_M
        elif self.dlg.cbDistUOM.currentIndex() == 1:  # KM
            dist_uom = UOM_KM
        elif self.dlg.cbDistUOM.currentIndex() == 2:  # NM
            dist_uom = UOM_NM
        elif self.dlg.cbDistUOM.currentIndex() == 3:  # feet
            dist_uom = UOM_FEET
        elif self.dlg.cbDistUOM.currentIndex() == 4:  # SM
            dist_uom = UOM_SM
        return dist_uom

    def create_mlyr(self, lyr_name):
        """ Create temporary 'memory' layer to store results of calculations
        :param lyr_name: string, layer name
        """
        tmp_lyr = QgsVectorLayer('Polygon?crs=epsg:4326', lyr_name, 'memory')
        prov = tmp_lyr.dataProvider()
        tmp_lyr.startEditing()
        prov.addAttributes([QgsField("POL_NAME", QVariant.String)])
        tmp_lyr.commitChanges()
        QgsMapLayerRegistry.instance().addMapLayers([tmp_lyr])

    def create_rectangle_buffer(self):
        if self.check_input():
            self.polygon_name = self.dlg.lePolygonName.text()
            buffer_m = to_meters(self.buffer, self.get_buffer_uom())
            dist_am = to_meters(self.dist_a, self.get_dist_uom())
            dist_bm = to_meters(self.dist_b, self.get_dist_uom())

            p1_lat, p1_lon = dist_azm_orth_offset2latlon(self.ref_lat, self.ref_lon, self.brng, dist_am, buffer_m, "RIGHT")
            p4_lat, p4_lon = dist_azm_orth_offset2latlon(self.ref_lat, self.ref_lon, self.brng, dist_am, buffer_m, "LEFT")

            p2_lat, p2_lon = dist_azm_orth_offset2latlon(self.ref_lat, self.ref_lon, self.brng, dist_bm, buffer_m, "RIGHT")
            p3_lat, p3_lon = dist_azm_orth_offset2latlon(self.ref_lat, self.ref_lon, self.brng, dist_bm, buffer_m, "LEFT")

            v1 = QgsPoint(p1_lon, p1_lat)
            v2 = QgsPoint(p2_lon, p2_lat)
            v3 = QgsPoint(p3_lon, p3_lat)
            v4 = QgsPoint(p4_lon, p4_lat)

            vertices = [v1, v2, v3, v4]

            feat = QgsFeature()
            feat.setGeometry(QgsGeometry.fromPolygon([vertices]))
            feat.setAttributes([self.polygon_name])

            layers = self.iface.legendInterface().layers()
            layer_list = []  # List of layers in current (opened) QGIS project

            for layer in layers:
                layer_list.append(layer.name())

            if self.mlyr_name not in layer_list:
                self.mlyr_name = "BUFFER_" + get_tmp_name()
                self.create_mlyr(self.mlyr_name)
                mlyr = QgsVectorLayer('Polygon?crs=epsg:4326', self.mlyr_name, 'memory')
                mlyr = self.iface.activeLayer()
                mlyr.startEditing()
                mprov = mlyr.dataProvider()
                mprov.addFeatures([feat])
                mlyr.commitChanges()
                mlyr.updateExtents()
                self.iface.mapCanvas().refresh()
            elif self.mlyr_name in layer_list:
                mlyr = QgsVectorLayer('Polygon?crs=epsg:4326', self.mlyr_name, 'memory')
                mlyr = self.iface.activeLayer()
                mlyr.startEditing()
                mprov = mlyr.dataProvider()
                mprov.addFeatures([feat])
                mlyr.commitChanges()
                mlyr.updateExtents()
                self.iface.mapCanvas().refresh()
        return

    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass
