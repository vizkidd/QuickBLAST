#include <string>
#include <vector>

#if defined(_WIN32) || defined(__MINGW32__) || defined(MINGW32) || defined(WIN32)
// #define _WIN32_WINNT _WIN32_WINNT_WIN7
#ifdef QBLIBRARY_EXPORTS
#define QBLIBRARY_API __declspec(dllexport)
#endif
#else
#define QBLIBRARY_API //__declspec(dllimport)
#endif

struct FastaSequenceData
{
    int rec_no = 1;
    std::string header;
    std::string seq;
};
struct BLASTHitData
{
    int rec_no = 1;
    std::vector<std::string_view> col_names();
    std::vector<std::string_view> col_values();
};