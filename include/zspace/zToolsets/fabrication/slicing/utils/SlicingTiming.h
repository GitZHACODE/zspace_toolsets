#ifndef ZSPACE_SLICING_TIMING_H
#define ZSPACE_SLICING_TIMING_H

#include <zspace/zToolsets/zToolsetsExport.h>
#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

namespace zSpace {
namespace zTsUtils {

// Inclusive elapsed times. Worker records sum elapsed time over calls and may
// exceed the containing parallel method's wall time; do not add nested rows.
struct SlicingTimingRecord {
  std::string name;
  std::uint64_t calls = 0;
  double totalMilliseconds = 0;
  double maxMilliseconds = 0;
  bool workerAggregate = false;
};

// Disabled by default. Enable/clear around one pipeline invocation, wait for its
// workers to finish, then read and disable. Recording is safe across workers.
ZSPACE_TOOLSETS bool slicingTimingsEnabled() noexcept;
ZSPACE_TOOLSETS void setSlicingTimingsEnabled(bool enabled) noexcept;
ZSPACE_TOOLSETS void clearSlicingTimings();
ZSPACE_TOOLSETS std::vector<SlicingTimingRecord> getSlicingTimings();

class ZSPACE_TOOLSETS ScopedSlicingTiming {
public:
  explicit ScopedSlicingTiming(const char *name, bool workerAggregate = false) noexcept;
  ~ScopedSlicingTiming() noexcept;
  void stop() noexcept;
  ScopedSlicingTiming(const ScopedSlicingTiming &) = delete;
  ScopedSlicingTiming &operator=(const ScopedSlicingTiming &) = delete;

private:
  const char *_name;
  bool _workerAggregate;
  std::uint64_t _generation;
  std::chrono::steady_clock::time_point _start;
};

} // namespace zTsUtils
} // namespace zSpace
#endif
