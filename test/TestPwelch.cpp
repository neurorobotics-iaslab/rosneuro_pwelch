#include "rosneuro_pwelch/Pwelch.hpp"
#include <gtest/gtest.h>

namespace rosneuro {

class TestPwelchSuite : public ::testing::Test {
    public:
        TestPwelchSuite() {}
        ~TestPwelchSuite() {}

        void SetUp() override {
            pwelch = new Pwelch<double>();
        }

        Pwelch<double>* pwelch;
        unsigned int wlength = 4;
        unsigned int novl = 0;
        unsigned int fs = 1000;
        unsigned int dolog = 0;
};

TEST_F(TestPwelchSuite, Constructor) {
    EXPECT_FALSE(pwelch->isset_);
    EXPECT_FALSE(pwelch->psd_configured_);
}

TEST_F(TestPwelchSuite, Configure) {
    EXPECT_TRUE(pwelch->configure(wlength, 0, novl, fs, dolog));
    EXPECT_EQ(pwelch->config.wlength, wlength);
    EXPECT_EQ(pwelch->config.novl, novl);
    EXPECT_EQ(pwelch->config.fs, fs);
    EXPECT_EQ(pwelch->config.dolog, dolog);
    EXPECT_EQ(pwelch->config.grid, pwelch->compute_grid(wlength, pwelch->config.nfft, fs));
    std::vector<int> wtypes = {1, 2, 3, 4};
    for (auto wtype : wtypes) {
        EXPECT_TRUE(pwelch->configure(wlength, wtype, novl, fs, dolog));
        EXPECT_EQ(pwelch->config.wtype, wtype);
        EXPECT_TRUE(pwelch->isset_);
        EXPECT_FALSE(pwelch->psd_configured_);
    }
}

TEST_F(TestPwelchSuite, Isset) {
    pwelch->isset_ = false;
    EXPECT_THROW(pwelch->isSet(), std::runtime_error);
}

TEST_F(TestPwelchSuite, Apply) {
    pwelch->psd_configured_ = false;

    DynamicMatrix<double> in(4, 2);
    in << 0.68, 0.82, -0.21, -0.60,
          0.56, -0.32, 0.59, 0.53;

    std::vector<DynamicMatrix<double>> matrix_list{
            (DynamicMatrix<double>(2, 2) << -3.10, -2.85, -4.49, -4.23).finished(),
            (DynamicMatrix<double>(2, 2) << -3.08, -2.93, -4.18, -3.95).finished(),
            (DynamicMatrix<double>(2, 2) << -3.08, -2.82, -4.47, -4.21).finished(),
            (DynamicMatrix<double>(2, 2) << -3.08, -2.82, -4.47, -4.21).finished()
    };

    for (int wtype = 1; wtype <= 4; ++wtype) {
        pwelch->configure(3, wtype, 2, 4, 1);
        EXPECT_TRUE(pwelch->apply(in).isApprox(matrix_list[wtype-1], 1e-2));
    }
}

}

int main(int argc, char **argv) {
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Fatal);
    ros::init(argc, argv, "test_pwelch");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}