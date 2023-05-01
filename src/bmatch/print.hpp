#ifndef PRINT_HPP
#define PRINT_HPP

#include <iostream>
#include <type_traits>

template <typename T, typename = void>
struct has_begin_end : std::false_type {};

template <typename T>
struct has_begin_end<T, decltype((void)std::begin(std::declval<T>()), (void)std::end(std::declval<T>()))> : std::true_type {};

template <typename T>
struct is_container : std::integral_constant<bool, has_begin_end<T>::value> {};

template <typename ...T>
struct is_container<std::basic_string<T...> > : std::false_type {};

template <typename T>
std::enable_if_t<is_container<T>::value, std::ostream&> operator<<(std::ostream& os, const T& container) {
    auto iter = std::begin(container);
    const auto end = std::end(container);

    os << '[';
    if (iter != end) {
        for ( ; ; ) {
            os << *iter;
            if (++iter == end) break;
            os << ", ";
        }
    }
    os << ']';

    return os;
}

template <typename ...T>
std::ostream& operator<<(std::ostream& os, const std::pair<T...>& p) {
    return os << "(" << p.first << ", " << p.second << ")";
}

static void print() {
    std::cout << std::endl;
}

template <typename T, typename ...Tail>
static void print(const T& t, Tail... tail) {
    std::cout << t << ' ';
    print(tail...);
}

#endif // PRINT_HPP
