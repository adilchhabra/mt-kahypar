include(FetchContent)
FetchContent_Populate(
  tbb
  URL https://github.com/uxlfoundation/oneTBB/releases/download/v2022.0.0/oneapi-tbb-2022.0.0-lin.tgz
  URL_HASH SHA256=1b669eb357dd90f3135f27e3c9a78683c6ecc74edf2799f7cb7df92a5423cb76
  SOURCE_DIR external_tools/tbb
)