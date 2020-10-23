// using standard exceptions
#include <iostream>
#include <exception>
using namespace std;

class myexception: public exception
{
  virtual const char* what() const throw()
  {
    return "My exception happened";
  }
} myex;

class yourexception: public exception
{
  virtual const char* what() const throw()
  {
    return "My exception tow  happened";
  }
} yourex;


int main () {
  try {
    throw yourex;
    throw myex;

  } 
catch (exception& e)
  {
    cout << e.what() << '\n';
    throw;
  }

  std::cout << "So what " << std::endl;
  return 0;
}
