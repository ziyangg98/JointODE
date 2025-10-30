// Stub implementation for CppAD functions not needed in R package
#include <string>
#include <stdexcept>

namespace CppAD {
namespace local {

// Provide a stub for temp_file() that is never actually called
// This is needed because some CppAD headers reference it, but we don't use
// the functionality that requires it (sparse Hessian computation to files)
std::string temp_file() {
    // This should never be called in our usage
    // If it is called, R will throw an error which is appropriate
    throw std::runtime_error("CppAD temp_file() not supported in R package");
    return "";
}

} // namespace local
} // namespace CppAD
