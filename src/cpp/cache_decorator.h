// From https://martin-ueding.de/posts/c-cache-decorator/.

#ifndef CACHE_DECORATOR_H
#define CACHE_DECORATOR_H

#include <functional>
#include <map>

template <typename R, typename... A>
class CacheDecorator {
  public:
    CacheDecorator(std::function<R(A...)> f) : f_(f) {}

    R operator()(A... a) {
        std::tuple<A...> key(a...);
        auto search = map_.find(key);
        if (search != map_.end()) {
            return search->second;
        }

        auto result = f_(a...);
        map_[key] = result;
        return result;
    }

  private:
    std::function<R(A...)> f_;
    std::map<std::tuple<A...>, R> map_;
};

#endif // CACHE_DECORATOR_H