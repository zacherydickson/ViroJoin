#ifndef SURVEYOR_STRUTILS_H
#define SURVEYOR_STRUTILS_H

#include <vector>
#include <string>
#include <sstream>
#include <type_traits>
#include <iterator>

//From Arafat Hasan: Answer to stack Overflow Question 14265581
std::vector<std::string> strsplit (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

//strjoin and to_strjoin written with 
//Templating help from gemini 3.5 Flash on May 29 2026
//  to ensure that only containers of strings or containers convertible to strings
//  can be used
template <typename It, 
          typename = typename std::enable_if<
              std::is_same<typename std::iterator_traits<It>::value_type, std::string>::value
          >::type>
std::string strjoin(It begin, It end, char delim) {
    std::string str = *begin;
    for (auto it = std::next(begin); it != end; ++it) {
        str += delim + *it;
    }
    return str;
}

//Functor serving as the default to_string
struct DefaultToString {
    template<typename T>
        auto operator()(const T& val) const -> decltype(std::to_string(val)) {
            return std::to_string(val);
        }
};

//USer can specify any functor/lambda they want as long as it produces a string
template <typename It, 
          typename F = DefaultToString,
          typename = typename std::enable_if<
              std::is_convertible<
                  decltype(std::declval<F>()(*std::declval<It>())),
                  std::string
              >::value
          >::type>
std::string to_strjoin(It begin, It end, char delim, F transformer = F()) {
    std::string str = transformer(*begin);
    for (auto it = std::next(begin); it != end; ++it) {
        str += delim + transformer(*it);
    }
    return str;
}

std::string to_upper(const std::string & str) {
    std::string out = str;
    for(auto & c : out) { c = (char)toupper(c); }
    return out;
}

#endif
