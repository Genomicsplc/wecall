package(default_visibility = ["//visibility:public"])

cc_library(
    name = "libgenomics.so",
    srcs = ["BUILD"] + glob([
        "src/**/*.cpp"
    ]),
    hdrs = ["BUILD"] + glob([
        "src/**/*.hpp",
    ]),
    linkshared = 1,
)

cc_binary(
    name = "weCall",
    srcs = ["src/weCall.cpp"],
    deps = [":libgenomics.so"],
)

filegroup(
    name = "srcs",
    srcs = ["BUILD"] + glob([
        "**/*.cpp",
        "**/*.hpp",
    ]),
)
