#include <zspace/zToolsets/fabrication/slicing/utils/SlicingTiming.h>
#include <algorithm>
#include <atomic>
#include <map>
#include <mutex>
#include <utility>

namespace zSpace {
namespace zTsUtils {
namespace {
std::atomic<bool> enabled{false};
std::atomic<std::uint64_t> generation{0};
std::mutex timingMutex;
std::map<std::pair<std::string, bool>, SlicingTimingRecord> timings;
}

bool slicingTimingsEnabled() noexcept { return enabled.load(std::memory_order_relaxed); }
void setSlicingTimingsEnabled(bool value) noexcept { enabled.store(value, std::memory_order_relaxed); }
void clearSlicingTimings() {
  std::lock_guard<std::mutex> lock(timingMutex);
  generation.fetch_add(1, std::memory_order_relaxed);
  timings.clear();
}
std::vector<SlicingTimingRecord> getSlicingTimings() {
  std::lock_guard<std::mutex> lock(timingMutex);
  std::vector<SlicingTimingRecord> result;
  result.reserve(timings.size());
  for (const auto &entry : timings) result.push_back(entry.second);
  return result;
}
ScopedSlicingTiming::ScopedSlicingTiming(const char *name, bool workerAggregate) noexcept
    : _name(slicingTimingsEnabled() ? name : nullptr),
      _workerAggregate(workerAggregate),
      _generation(generation.load(std::memory_order_relaxed)) {
  if (_name) _start = std::chrono::steady_clock::now();
}
ScopedSlicingTiming::~ScopedSlicingTiming() noexcept { stop(); }
void ScopedSlicingTiming::stop() noexcept {
  if (!_name) return;
  const char *name = _name;
  _name = nullptr;
  const double elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - _start).count();
  // Instrumentation must never change success/failure of a geometry operation.
  try {
    std::lock_guard<std::mutex> lock(timingMutex);
    if (!slicingTimingsEnabled() || _generation != generation.load(std::memory_order_relaxed)) return;
    auto &record = timings[{name, _workerAggregate}];
    record.name = name;
    record.workerAggregate = _workerAggregate;
    ++record.calls;
    record.totalMilliseconds += elapsed;
    record.maxMilliseconds = (std::max)(record.maxMilliseconds, elapsed);
  } catch (...) {}
}

} // namespace zTsUtils
} // namespace zSpace
