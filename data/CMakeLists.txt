# Tracking geometry
INSTALL(FILES MuColl_v1.root MuSIC_v1.root MuSIC_v2.root
  DESTINATION ${CMAKE_INSTALL_DATADIR}/MarlinACTS/data)

# Geometry description
INSTALL(FILES MuColl_v1.json MuSIC_v1.json MuSIC_v2.json
  DESTINATION ${CMAKE_INSTALL_DATADIR}/MarlinACTS/data)

# Material map
INSTALL(FILES material-maps.json
  DESTINATION ${CMAKE_INSTALL_DATADIR}/MarlinACTS/data)
