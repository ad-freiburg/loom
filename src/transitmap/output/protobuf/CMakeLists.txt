if (Protobuf_FOUND AND Protobuf_PROTOC_EXECUTABLE)
include_directories(${Protobuf_INCLUDE_DIR})

protobuf_generate_cpp(PROTO_SRC PROTO_HEADER vector_tile.proto)
add_library(vt-proto ${PROTO_HEADER} ${PROTO_SRC})

target_link_libraries(vt-proto PUBLIC protobuf::libprotobuf-lite)
endif()
