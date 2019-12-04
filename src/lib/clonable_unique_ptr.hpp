#pragma once

template <typename T>
struct ClonableUniquePtr {
    std::unique_ptr<T> ptr;

    ClonableUniquePtr() = default;
    ClonableUniquePtr(ClonableUniquePtr &&) = default;
    ClonableUniquePtr(std::unique_ptr<T> &&other) : ptr{std::move(other)} {};
    ClonableUniquePtr(ClonableUniquePtr const &other) : ptr{other.ptr->clone()} {};
    T *operator->() const { return ptr; }
};