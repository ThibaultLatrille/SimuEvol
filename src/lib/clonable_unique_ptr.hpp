#pragma once

template <typename T>
struct ClonableUniquePtr {
  private:
    std::unique_ptr<T> ptr;

  public:
    T *get() const { return ptr.get(); }
    ClonableUniquePtr() = default;
    ClonableUniquePtr(ClonableUniquePtr &&) = default;
    ClonableUniquePtr(std::unique_ptr<T> &&other) : ptr{std::move(other)} {};
    ClonableUniquePtr(ClonableUniquePtr const &other) : ptr{other.ptr->clone()} {};
    ClonableUniquePtr &operator=(ClonableUniquePtr const &other) {
        ptr.reset(other.ptr->clone());
        return *this;
    }
    ClonableUniquePtr &operator=(ClonableUniquePtr &&other) = default;
    T *operator->() const { return ptr.get(); }
    T &operator*() const { return *ptr; }
};