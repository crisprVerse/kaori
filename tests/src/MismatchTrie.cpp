#include <gtest/gtest.h>
#include "kaori/MismatchTrie.hpp"
#include <string>

TEST(MismatchTrie, Basic) {
    std::vector<std::string> things { "ACGT", "AAAA", "ACAA", "AGTT" };

    std::vector<const char*> ptrs;
    for (const auto& t : things) {
        ptrs.push_back(t.c_str());
    }

    kaori::MismatchTrie stuff(ptrs, 4);
    {
        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 0);
    }

    {
        auto res = stuff.search("AAAT", 1);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("CCAG", 2);
        EXPECT_EQ(res.first, 2);
        EXPECT_EQ(res.second, 2);
    }

    {
        auto res = stuff.search("AGTT", 0);
        EXPECT_EQ(res.first, 3);
        EXPECT_EQ(res.second, 0);
    }
}

TEST(MismatchTrie, ReverseComplement) {
    std::vector<std::string> things { "ACGT", "TTTT", "TTGT", "AACT" };

    std::vector<const char*> ptrs;
    for (const auto& t : things) {
        ptrs.push_back(t.c_str());
    }

    kaori::MismatchTrie stuff(ptrs, 4, true);
    {
        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 0);
    }

    {
        auto res = stuff.search("AAAT", 1);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("CCAG", 2);
        EXPECT_EQ(res.first, 2);
        EXPECT_EQ(res.second, 2);
    }

    {
        auto res = stuff.search("AGTT", 0);
        EXPECT_EQ(res.first, 3);
        EXPECT_EQ(res.second, 0);
    }
}

TEST(MismatchTrie, MoreMismatches) {
    std::vector<std::string> things { "ACGTACGTACGT", "TTTGGGCCCAAA" };

    std::vector<const char*> ptrs;
    for (const auto& t : things) {
        ptrs.push_back(t.c_str());
    }

    kaori::MismatchTrie stuff(ptrs, 12);
    {
        auto res = stuff.search("ACGTACGTCCGT", 2);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("TCGTACGTCCGT", 2);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 2);
    }

    {
        auto res = stuff.search("TTTGGGGCCAAA", 2);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("TTGGGGGCCAAA", 2);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 2);
    }
}

TEST(MismatchTrie, CappedMismatch) {
    std::vector<std::string> things { "ACGT", "AAAA", "ACAA", "AGTT" };

    std::vector<const char*> ptrs;
    for (const auto& t : things) {
        ptrs.push_back(t.c_str());
    }

    kaori::MismatchTrie stuff(ptrs, 4);
    {
        auto res = stuff.search("AAAT", 0);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("CCAG", 1);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 2);
    }

    {
        auto res = stuff.search("CCAG", 1);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 2);
    }
}

TEST(MismatchTrie, Ambiguous) {
    std::vector<std::string> things { "AAAAGAAAA", "AAAACAAAA", "AAAAAAAAG", "AAAAAAAAC" };

    std::vector<const char*> ptrs;
    for (const auto& t : things) {
        ptrs.push_back(t.c_str());
    }

    kaori::MismatchTrie stuff(ptrs, 9);
    {
        auto res = stuff.search("AAAATAAAA", 1);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("AAAAAAAAT", 1);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 1);
    }
}

