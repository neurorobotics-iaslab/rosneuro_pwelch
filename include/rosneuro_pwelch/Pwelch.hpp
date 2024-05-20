#ifndef ROSNEURO_PWELCH_PWELCH_HPP
#define ROSNEURO_PWELCH_PWELCH_HPP

#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <fftw3.h>
#include <ros/ros.h>
#include <rosneuro_msgs/NeuroFrame.h>
#include <gtest/gtest_prod.h>
#include "rosneuro_windows_hann/Hann.hpp"
#include "rosneuro_windows_hamming/Hamming.hpp"
#include "rosneuro_windows/Window.hpp"
#include "rosneuro_windows_blackman/Blackman.hpp"
#include "rosneuro_windows_flattop/Flattop.hpp"

namespace rosneuro {
    typedef struct {
        unsigned int wlength;
        unsigned int wtype;
        unsigned int novl;
        unsigned int nfft;
        unsigned int fs;
        unsigned int dolog;
        std::vector<uint32_t> 	grid;
    } pwelch_config;

    template<typename T> using DynamicMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    template <typename T>
    class Pwelch {
        public:
            Pwelch(void);
            ~Pwelch(void);

            bool configure(unsigned int wlength, unsigned int wtype, unsigned novl, unsigned int fs, unsigned int dolog);
            DynamicMatrix<T> apply(const DynamicMatrix<T>& in);
            bool isSet(void);
            pwelch_config config;
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        private:
            unsigned int compute_nfft(unsigned int wlength);
            unsigned int compute_nsegments(unsigned int wlength, unsigned int novl, unsigned int nsamples);
            void compute(const Eigen::Ref<const Eigen::VectorXd>& in, Eigen::Ref<Eigen::VectorXd> out);
            std::vector<uint32_t> compute_grid(unsigned int wlength, unsigned int nfft, unsigned int fs);
            void setWindow(void);

            rosneuro_msgs::NeuroFrame neuromsg_;
            Window<T> * window_;
            Eigen::VectorXd  wsig_;
            Eigen::VectorXd  wpxx_;
            fftw_plan 	plan_;
            DynamicMatrix<T>  psd_;
            bool 		isset_;
            bool        psd_configured_;

            FRIEND_TEST(TestPwelchSuite, Constructor);
            FRIEND_TEST(TestPwelchSuite, Configure);
            FRIEND_TEST(TestPwelchSuite, Isset);
            FRIEND_TEST(TestPwelchSuite, Apply);
    };

    template <typename T>
    Pwelch<T>::Pwelch(void){
        this->isset_ = false;
        this->psd_configured_ = false;
    }

    template <typename T>
    Pwelch<T>::~Pwelch(void) {
        fftw_free(this->plan_);
        free(this->window_);
    }

    template <typename T>
    bool Pwelch<T>::configure(unsigned int wlength, unsigned int wtype, unsigned novl, unsigned int fs, unsigned int dolog) {
        this->config.wlength	= wlength;
        this->config.wtype 	= wtype;
        this->config.novl	= novl;
        this->config.nfft 	= compute_nfft(config.wlength);
        this->config.fs 	= fs;
        this->config.dolog 	= dolog;
        this->config.grid 	= compute_grid(config.wlength, this->config.nfft, config.fs);

        this->setWindow();

        this->wsig_ 	= Eigen::VectorXd::Zero(this->config.wlength);
        this->wpxx_ 	= Eigen::VectorXd::Zero(this->config.wlength);
        this->plan_ = fftw_plan_r2r_1d(this->config.wlength, this->wsig_.data(), this->wpxx_.data(), FFTW_R2HC, FFTW_PATIENT);
        if(this->plan_ == NULL){
            return false;
        }

        this->psd_ = Eigen::MatrixXd::Zero(this->config.nfft, 1);
        this->isset_ = true;
        this->psd_configured_ = false;
        return this->isset_;
    }

    template <typename T>
    void Pwelch<T>::setWindow(){
        switch (this->config.wtype){
            case 1:
                this->window_ = new Hann<double>();
                break;
            case 2:
                this->window_ = new Hamming<double>();
                break;
            case 3:
                this->window_ = new Blackman<double>();
                break;
            case 4:
                this->window_ = new Flattop<double>();
                break;
            default:
                break;
        }
    }

    template <typename T>
    bool Pwelch<T>::isSet(void){
        if(!this->isset_) {
            ROS_ERROR("The pwelch is not configured");
            throw std::runtime_error("The pwelch is not configured");
        }
        return this->isset_;
    }

    template <typename T>
    DynamicMatrix<T> Pwelch<T>::apply(const DynamicMatrix<T>& in){
        unsigned int nchannels = in.cols();

        if(!this->psd_configured_){
            this->psd_.resize(Eigen::NoChange, in.cols());
            this->psd_configured_ = true;
            ROS_WARN("[Pwelch] First apply: the psd size is set");
        }

        for (unsigned int i = 0; i < nchannels; i++) {
            compute(in.col(i), this->psd_.col(i));
        }

        if (this->config.dolog) {
            this->psd_ = this->psd_.array().log();
        }

        return this->psd_;
    }

    template <typename T>
    void Pwelch<T>::compute(const Eigen::Ref<const Eigen::VectorXd>& in, Eigen::Ref<Eigen::VectorXd> out) {
        unsigned int nsamples = in.size();
        unsigned int wlength = this->config.wlength;
        unsigned int novl = this->config.novl;
        unsigned int nfft = this->config.nfft;
        unsigned int fs = this->config.fs;

        int nsegments = compute_nsegments(wlength, novl, nsamples);
        Eigen::VectorXd power_spectral_density	= Eigen::VectorXd::Zero(nfft);

        for (unsigned int segId = 0; segId < nsegments; ++segId) {
            this->wsig_ = this->window_->apply(in.segment(segId*(wlength - novl), wlength));
            fftw_execute_r2r(this->plan_, this->wsig_.data(), this->wpxx_.data());

            // out spans from 0 to wLength/2 (NFFT = wLength/2 + 1)
            power_spectral_density(0) += pow(this->wpxx_(0), 2);
            power_spectral_density(wlength/2) += pow(this->wpxx_(wlength/2), 2);

            //	out(k) 	 = wpxx(k)^2 + wpxx(n-k)^2, k = 1 : (n/2 - 1) [Real and imagery part]
            power_spectral_density.segment(1, wlength/2 - 1).array() += this->wpxx_.segment(1, wlength/2 - 1).array().pow(2)
                                       + this->wpxx_.segment(wlength/2 + 1, wlength/2 - 1).reverse().array().pow(2);

        }

        double window_norm = this->window_->GetWindowNorm();
        double power_spectral_density_factor = (nsegments * window_norm * fs * wlength)/2.0;

        power_spectral_density(0) /= (2.0 * power_spectral_density_factor);
        power_spectral_density(wlength / 2) /= (2.0 * power_spectral_density_factor);
        power_spectral_density.segment(1, wlength / 2 - 1).array() /= power_spectral_density_factor;

        out = power_spectral_density;
    }

    template <typename T>
    unsigned int Pwelch<T>::compute_nsegments(unsigned int wlength, unsigned int novl, unsigned int nsamples) {

        if (wlength == novl)
            novl = wlength - 1;

        return floor((nsamples - novl) / (wlength - novl));
    }

    template <typename T>
    unsigned int Pwelch<T>::compute_nfft(unsigned int wlength) {
        return floor(wlength/2) + 1;
    }

    template <typename T>
    std::vector<uint32_t> Pwelch<T>::compute_grid(unsigned int wlength, unsigned int nfft, unsigned int fs) {
        std::vector<uint32_t> vgrid;
        Eigen::VectorXi grid = Eigen::VectorXi::LinSpaced(nfft, 0, nfft - 1) * (int)std::floor(fs/wlength);
        for (unsigned int i=0; i< grid.size(); i++) {
            vgrid.push_back(grid(i));
        }
        return vgrid;
    }
}

#endif