#include <gtest/gtest.h>
#include "kaori/handlers/DualBarcodes.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class DualBarcodesTest : public testing::Test {
protected:
    DualBarcodesTest() : 
        constant1("AAAA----CGGC"),
        constant2("AGCT------TTTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::vector<std::vector<const char*> > make_pointers() const {
        return std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) };
    }

    std::pair<const char*, const char*> bounds(const std::string& s) {
        return std::make_pair(s.c_str(), s.c_str() + s.size());
    }

    std::string constant1, constant2;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;
};

TEST_F(DualBarcodesTest, BasicFirst) {
    kaori::DualBarcodes<128> stuff(
        constant1.c_str(), constant1.size(), false, to_pointers(variables1), 0,
        constant2.c_str(), constant2.size(), false, to_pointers(variables2), 0
    );
}
