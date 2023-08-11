#include <type_traits>
#include <iostream>

int main() {
    std::cout << "hello\n";
    std::cout << std::is_constructible<int>::value << "\n";
    std::cout << "hello\n";
    std::cout << std::is_default_constructible<int>::value << "\n";
    std::cout << "hello\n";
    std::cout << std::is_trivially_copyable<int>::value << "\n"; 
    std::cout << "hello\n";
   if (__cplusplus == 201703L)
      std::cout << "C++17" << std::endl;
   else if (__cplusplus == 201402L)
      std::cout << "C++14" << std::endl;
   else if (__cplusplus == 201103L)
      std::cout << "C++11" << std::endl;
   else if (__cplusplus == 199711L)
      std::cout << "C++98" << std::endl;
   else
      std::cout << "pre-standard C++" << std::endl;

}
