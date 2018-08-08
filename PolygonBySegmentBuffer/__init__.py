# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PolygonBySegmentBuffer
                                 A QGIS plugin
 Creates polygon which is defined as buffer around segment
                             -------------------
        begin                : 2018-08-08
        copyright            : (C) 2018 by Paweł Strzelewicz
        email                : pawel.strzelewicz83@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load PolygonBySegmentBuffer class from file PolygonBySegmentBuffer.

    :param iface: A QGIS interface instance.
    :type iface: QgisInterface
    """
    #
    from .polygon_by_segmentbuffer import PolygonBySegmentBuffer
    return PolygonBySegmentBuffer(iface)
