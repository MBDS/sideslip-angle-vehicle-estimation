
/* Example invokation:

./fg-estimate-slip-angle \
  --race-dataset-directory ../../datasets/20140222_01_01_03_250lm  \
  --sigma-factor-ay 7 \
  --sigma-obs-rate 1e-2 \
  --sigma-factor-rate 0.9e-2 \
  --sigma-factor-beta 4e-3

*/

#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/NonlinearEquality.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/factorTesting.h>  // for unit tests
#include <gtsam/slam/PriorFactor.h>
#include <mrpt/3rdparty/tclap/CmdLine.h>
#include <mrpt/core/bits_math.h>
#include <mrpt/math/CMatrixDynamic.h>
#include <mrpt/math/CVectorDynamic.h>
#include <mrpt/system/CTimeLogger.h>

#include <fstream>
#include <iostream>

// Declare the supported command line switches ===========
TCLAP::CmdLine cmd("fg-estimate-slip-angle", ' ');

TCLAP::ValueArg<std::string> arg_race_dir(
    "", "race-dataset-directory", "Directory with race dataset", true,
    "race_xxx", "Race dataset directory", cmd);

TCLAP::ValueArg<double> arg_sigma_factor_beta(
    "", "sigma-factor-beta", "Noise model", false, 4e-3, "Sigma", cmd);
TCLAP::ValueArg<double> arg_sigma_factor_rate(
    "", "sigma-factor-rate", "Noise model", false, 0.009, "Sigma", cmd);
TCLAP::ValueArg<double> arg_sigma_obs_rate(
    "", "sigma-obs-rate", "Noise model", false, 1e-2, "Sigma", cmd);
TCLAP::ValueArg<double> arg_sigma_factor_ay(
    "", "sigma-factor-ay", "Noise model", false, 7.0, "Sigma", cmd);

TCLAP::ValueArg<int> arg_fixed_lag_window(
    "", "fixed-lag-window",
    "If ==0, use batch optimizer. If >=1, use fixed-lag smoother", false, 0,
    "WindowLength", cmd);

TCLAP::ValueArg<std::string> arg_output_prefix(
    "", "output-prefix", "Output files prefix", false, "out_", "PREFIX", cmd);

struct CarParameters
{
    const double lf    = 1.33;  // m
    const double lr    = 1.07;  // m
    const double m     = 982;  // kg
    const double a     = lf + 0.841;
    const double b     = lr + 0.849;  // m
    const double h     = 1.115;  // m
    const double track = 1.35;  // m  (che poi sarebbe tf e/o tr)
    const double Len   = 4.09;  // m
    const double Jz    = m / 12 * (1.7 * 1.7 + 4.09 * 4.09);  // kgm^2
    const double L     = lr + lf;  // m
    const double tf    = 1.35;  // m
    const double tr    = tf;  // m

    const double Cf = 0.7e5;
    const double Cr = 1.2e5;
};

struct CarDataset
{
    double                    dt;
    mrpt::math::CMatrixDouble t;
    mrpt::math::CMatrixDouble ax;
    mrpt::math::CMatrixDouble ay;
    mrpt::math::CMatrixDouble yawRate;  // [rad/s]
    mrpt::math::CMatrixDouble vx, vy;
    mrpt::math::CMatrixDouble delta;

    mrpt::math::CMatrixDouble beta_true;  // ground truth
};

CarDataset load_car_dataset(const std::string& directory)
{
    using namespace std::string_literals;

    CarDataset d;
    d.dt = 0.01;  // 100 Hz
    d.t.loadFromTextFile(directory + "/t.txt"s);
    d.ay.loadFromTextFile(directory + "/ayCG.txt"s);
    d.ax.loadFromTextFile(directory + "/axCG.txt"s);
    d.yawRate.loadFromTextFile(directory + "/yawRate.txt"s);
    d.vx.loadFromTextFile(directory + "/vx.txt"s);
    d.vy.loadFromTextFile(directory + "/vy.txt"s);
    d.delta.loadFromTextFile(directory + "/delta.txt"s);
    d.beta_true.loadFromTextFile(directory + "/beta_true.txt"s);

    ASSERT_EQUAL_(d.t.size(), d.ay.size());
    ASSERT_EQUAL_(d.t.size(), d.ax.size());
    ASSERT_EQUAL_(d.t.size(), d.yawRate.size());
    ASSERT_EQUAL_(d.t.size(), d.vx.size());
    ASSERT_EQUAL_(d.t.size(), d.vy.size());
    ASSERT_EQUAL_(d.t.size(), d.delta.size());
    ASSERT_EQUAL_(d.t.size(), d.beta_true.size());

    // Rfl = 0.2976;  // m
    // Rfr = 0.2976;  // m
    // Rrl = 0.3199;  // m
    // Rrr = 0.3199;  // m

    return d;
}

// ========= Factor definitions   =============

/// ==== FactorBeta ====
class FactorBeta
    : public gtsam::NoiseModelFactor3<
          double /* beta_{k-1} */, double /* beta_{k} */, double /* r_{k-1} */>
{
   private:
    using This = FactorBeta;
    using Base = gtsam::NoiseModelFactor3<
        double /* beta_{k-1} */, double /* beta_{k} */, double /* r_{k-1} */>;

    CarParameters carParameters_;
    double        dt_;  //!< Timestep [s]
    double        u_;  //!< Longitudinal velocity [m/s]
    double        delta_;  //!< Steering angle

   public:
    // shorthand for a smart pointer to a factor
    using shared_ptr = boost::shared_ptr<This>;

    /** default constructor - only use for serialization */
    FactorBeta() = default;

    /** Constructor */
    FactorBeta(
        const CarParameters& carParameters, const double dt, const double u,
        const double delta, const gtsam::SharedNoiseModel& noiseModel,
        gtsam::Key key_beta_k_1, gtsam::Key key_beta_k, gtsam::Key key_r_1)
        : Base(noiseModel, key_beta_k_1, key_beta_k, key_r_1),
          carParameters_(carParameters),
          dt_(dt),
          u_(u),
          delta_(delta)
    {
    }

    virtual ~FactorBeta() override {}

    /// @return a deep copy of this factor
    virtual gtsam::NonlinearFactor::shared_ptr clone() const override
    {
        return boost::static_pointer_cast<gtsam::NonlinearFactor>(
            gtsam::NonlinearFactor::shared_ptr(new This(*this)));
    }

    /** implement functions needed for Testable */

    /** print */
    virtual void print(
        const std::string& s, const gtsam::KeyFormatter& keyFormatter =
                                  gtsam::DefaultKeyFormatter) const override
    {
        std::cout << s << "FactorBeta(" << keyFormatter(key1()) << ","
                  << keyFormatter(key2()) << "," << keyFormatter(key3()) << ")"
                  << " dt: " << dt_ << " u: " << u_ << " delta: " << delta_
                  << "\n";

        noiseModel_->print("  noise model: ");
    }

    /** equals */
    virtual bool equals(
        const gtsam::NonlinearFactor& expected,
        double                        tol = 1e-9) const override
    {
        const This* e = dynamic_cast<const This*>(&expected);
        return e != nullptr && Base::equals(*e, tol);
    }

    gtsam::Vector evaluateError(
        const double& beta_k_1, const double& beta_k, const double& r_k_1,
        boost::optional<gtsam::Matrix&> J_beta_k_1,
        boost::optional<gtsam::Matrix&> J_beta_k,
        boost::optional<gtsam::Matrix&> J_r_k_1) const override
    {
        const auto& cp = carParameters_;

        // error:
        gtsam::Vector err(1);
        err[0] =
            beta_k - beta_k_1 -
            dt_ *
                (  //
                    -((cp.Cf + cp.Cr) / (cp.m * u_)) * beta_k_1  //
                    - ((cp.Cf * cp.lf - cp.Cr * cp.lr) / (cp.m * u_ * u_) + 1) *
                          r_k_1  //
                    + (cp.Cf * delta_ / (cp.m * u_)));

        MRPT_CHECK_NORMAL_NUMBER(err[0]);
        // error Jacobians:
        if (J_beta_k_1)
        {
            auto& J = J_beta_k_1.get();
            J.resize(1, 1);
            J(0, 0) = -(1.0 - dt_ * ((cp.Cf + cp.Cr) / (cp.m * u_)));
        }
        if (J_beta_k)
        {
            auto& J = J_beta_k.get();
            J.resize(1, 1);
            J(0, 0) = 1.0;
        }
        if (J_r_k_1)
        {
            auto& J = J_r_k_1.get();
            J.resize(1, 1);
            J(0, 0) =
                dt_ * ((cp.Cf * cp.lf - cp.Cr * cp.lr) / (cp.m * u_ * u_) + 1);
        }
        return err;
    }

   private:
    /** Serialization function */
    friend class boost::serialization::access;
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int /*version*/)
    {
        ar& boost::serialization::make_nvp(
            "FactorBeta", boost::serialization::base_object<Base>(*this));
    }
};
/// ==== FactorRate ====
class FactorRate
    : public gtsam::NoiseModelFactor3<
          double /* r_{k-1} */, double /* r_{k} */, double /* beta_{k-1} */>
{
   private:
    using This = FactorRate;
    using Base = gtsam::NoiseModelFactor3<
        double /* r_{k-1} */, double /* r_{k} */, double /* beta_{k-1} */>;

    CarParameters carParameters_;
    double        dt_;  //!< Timestep [s]
    double        u_;  //!< Longitudinal velocity [m/s]
    double        delta_;  //!< Steering angle

   public:
    // shorthand for a smart pointer to a factor
    using shared_ptr = boost::shared_ptr<This>;

    /** default constructor - only use for serialization */
    FactorRate() = default;

    /** Constructor */
    FactorRate(
        const CarParameters& carParameters, const double dt, const double u,
        const double delta, const gtsam::SharedNoiseModel& noiseModel,
        gtsam::Key key_r_k_1, gtsam::Key key_r_k, gtsam::Key key_beta_1)
        : Base(noiseModel, key_r_k_1, key_r_k, key_beta_1),
          carParameters_(carParameters),
          dt_(dt),
          u_(u),
          delta_(delta)
    {
    }

    virtual ~FactorRate() override {}

    /// @return a deep copy of this factor
    virtual gtsam::NonlinearFactor::shared_ptr clone() const override
    {
        return boost::static_pointer_cast<gtsam::NonlinearFactor>(
            gtsam::NonlinearFactor::shared_ptr(new This(*this)));
    }

    /** implement functions needed for Testable */

    /** print */
    virtual void print(
        const std::string& s, const gtsam::KeyFormatter& keyFormatter =
                                  gtsam::DefaultKeyFormatter) const override
    {
        std::cout << s << "FactorRate(" << keyFormatter(key1()) << ","
                  << keyFormatter(key2()) << "," << keyFormatter(key3()) << ")"
                  << " dt: " << dt_ << " u: " << u_ << " delta: " << delta_
                  << "\n";
        noiseModel_->print("  noise model: ");
    }

    /** equals */
    virtual bool equals(
        const gtsam::NonlinearFactor& expected,
        double                        tol = 1e-9) const override
    {
        const This* e = dynamic_cast<const This*>(&expected);
        return e != nullptr && Base::equals(*e, tol);
    }

    gtsam::Vector evaluateError(
        const double& r_k_1, const double& r_k, const double& beta_k_1,
        boost::optional<gtsam::Matrix&> J_r_k_1,
        boost::optional<gtsam::Matrix&> J_r_k,
        boost::optional<gtsam::Matrix&> J_beta_k_1) const override
    {
        const auto& cp = carParameters_;

        // error:

        const double lf = cp.lf, lr = cp.lr;
        const double lf2 = lf * lf, lr2 = lr * lr;

        gtsam::Vector err(1);
        err[0] =
            r_k - r_k_1 -
            dt_ * (  //
                      -((cp.Cf * lf - cp.Cr * lr) / cp.Jz) * beta_k_1  //
                      - ((cp.Cf * lf2 + cp.Cr * lr2) / (cp.Jz * u_)) * r_k_1  //
                      + (cp.Cf * lf * delta_) / cp.Jz);

        MRPT_CHECK_NORMAL_NUMBER(err[0]);
        // error Jacobians:
        if (J_r_k_1)
        {
            auto& J = J_r_k_1.get();
            J.resize(1, 1);
            J(0, 0) = -1.0 + dt_ * ((cp.Cf * lf2 + cp.Cr * lr2) / (cp.Jz * u_));
        }
        if (J_r_k)
        {
            auto& J = J_r_k.get();
            J.resize(1, 1);
            J(0, 0) = 1.0;
        }
        if (J_beta_k_1)
        {
            auto& J = J_beta_k_1.get();
            J.resize(1, 1);
            J(0, 0) = dt_ * ((cp.Cf * lf - cp.Cr * lr) / cp.Jz);
        }
        return err;
    }

   private:
    /** Serialization function */
    friend class boost::serialization::access;
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int /*version*/)
    {
        ar& boost::serialization::make_nvp(
            "FactorRate", boost::serialization::base_object<Base>(*this));
    }
};

/// ==== FactorLateralAcceleration ====
class FactorLateralAcceleration
    : public gtsam::NoiseModelFactor2<double /* beta_{k} */, double /* r_{k} */>
{
   private:
    using This = FactorLateralAcceleration;
    using Base =
        gtsam::NoiseModelFactor2<double /* beta_{k} */, double /* r_{k} */>;

    CarParameters carParameters_;
    double        u_;  //!< Longitudinal velocity [m/s]
    double        delta_;  //!< Steering angle
    double        ay_;  //!< Sensed lateral acceleration [m/s2]

   public:
    // shorthand for a smart pointer to a factor
    using shared_ptr = boost::shared_ptr<This>;

    /** default constructor - only use for serialization */
    FactorLateralAcceleration() = default;

    /** Constructor */
    FactorLateralAcceleration(
        const CarParameters& carParameters, const double u, const double delta,
        const double ay, const gtsam::SharedNoiseModel& noiseModel,
        gtsam::Key key_beta_k, gtsam::Key key_r_k)
        : Base(noiseModel, key_beta_k, key_r_k),
          carParameters_(carParameters),
          u_(u),
          delta_(delta),
          ay_(ay)
    {
    }

    virtual ~FactorLateralAcceleration() override {}

    /// @return a deep copy of this factor
    virtual gtsam::NonlinearFactor::shared_ptr clone() const override
    {
        return boost::static_pointer_cast<gtsam::NonlinearFactor>(
            gtsam::NonlinearFactor::shared_ptr(new This(*this)));
    }

    /** implement functions needed for Testable */

    /** print */
    virtual void print(
        const std::string& s, const gtsam::KeyFormatter& keyFormatter =
                                  gtsam::DefaultKeyFormatter) const override
    {
        std::cout << s << "FactorLateralAcceleration(" << keyFormatter(key1())
                  << "," << keyFormatter(key2()) << ")"
                  << " ay: " << ay_ << " u: " << u_ << " delta: " << delta_
                  << "\n";
        noiseModel_->print("  noise model: ");
    }

    /** equals */
    virtual bool equals(
        const gtsam::NonlinearFactor& expected,
        double                        tol = 1e-9) const override
    {
        const This* e = dynamic_cast<const This*>(&expected);
        return e != nullptr && Base::equals(*e, tol);
    }

    gtsam::Vector evaluateError(
        const double& beta_k, const double& r_k,
        boost::optional<gtsam::Matrix&> J_beta_k,
        boost::optional<gtsam::Matrix&> J_r_k) const override
    {
        const auto& cp = carParameters_;

        // error:
        const double lf = cp.lf, lr = cp.lr;

        gtsam::Vector err(1);
        err[0] = ay_ + ((cp.Cf + cp.Cr) / cp.m) * beta_k  //
                 + ((cp.Cf * lf - cp.Cr * lr) / (cp.m * u_)) * r_k  //
                 - (cp.Cf * delta_) / cp.m;

        MRPT_CHECK_NORMAL_NUMBER(err[0]);
        // error Jacobians:
        if (J_beta_k)
        {
            auto& J = J_beta_k.get();
            J.resize(1, 1);
            J(0, 0) = ((cp.Cf + cp.Cr) / cp.m);
        }
        if (J_r_k)
        {
            auto& J = J_r_k.get();
            J.resize(1, 1);
            J(0, 0) = ((cp.Cf * lf - cp.Cr * lr) / (cp.m * u_));
        }
        return err;
    }

   private:
    /** Serialization function */
    friend class boost::serialization::access;
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int /*version*/)
    {
        ar& boost::serialization::make_nvp(
            "FactorLateralAcceleration",
            boost::serialization::base_object<Base>(*this));
    }
};
// ====== End of factor definitions ===========

void slip_angle_smoother()
{
    mrpt::system::CTimeLogger profiler(true, "slip_angle_smoother");

    const auto sBeta    = gtsam::SymbolGenerator('b');
    const auto sYawRate = gtsam::SymbolGenerator('r');

    const auto carDataset = load_car_dataset(arg_race_dir.getValue());

    const CarParameters carParams;

    const size_t N = carDataset.t.size();

    std::cout << "Loaded race dataset '" << arg_race_dir.getValue() << "' with "
              << N << " entries.\n";

    mrpt::system::CTimeLoggerEntry tle1(profiler, "1.build_fg");

    const double dt = carDataset.dt;

    auto noise_prior_beta    = gtsam::noiseModel::Isotropic::Sigma(1, 1e2);
    auto noise_prior_yawrate = gtsam::noiseModel::Isotropic::Sigma(1, 1e2);

    auto noise_prior_window = gtsam::noiseModel::Isotropic::Sigma(1, 1.0);

    auto noise_factor_beta = gtsam::noiseModel::Isotropic::Sigma(
        1, arg_sigma_factor_beta.getValue());

    auto noise_factor_rate = gtsam::noiseModel::Isotropic::Sigma(
        1, arg_sigma_factor_rate.getValue());

    auto noise_observation_rate =
        gtsam::noiseModel::Isotropic::Sigma(1, arg_sigma_obs_rate.getValue());

    auto noise_factor_ay =
        gtsam::noiseModel::Isotropic::Sigma(1, arg_sigma_factor_ay.getValue());

    // Create Prior factors:
    double last_beta    = .0;
    double last_yawrate = .0;

    double t = .0;  // time

    gtsam::NonlinearFactorGraph wholeFG;
    gtsam::Values               wholeValues;

    wholeFG.emplace_shared<gtsam::PriorFactor<double>>(
        sBeta(0), last_beta, noise_prior_beta);

    wholeFG.emplace_shared<gtsam::PriorFactor<double>>(
        sYawRate(0), last_yawrate, noise_prior_yawrate);

    wholeValues.insert(sBeta(0), last_beta);
    wholeValues.insert(sYawRate(0), last_yawrate);

    // Save states to files:
    mrpt::math::CMatrixDouble all_beta(N, 2), all_yawrate(N, 2),
        all_obs_r(N, 2);

    auto lambda_Values_toAll = [&](const gtsam::Values& values) {
        for (const auto& kv : values)
        {
            gtsam::Symbol s(kv.key);
            const auto    step = s.index();
            const double  val  = values.at<double>(s);

            switch (s.chr())
            {
                case 'b':
                    all_beta(step, 1) = val;
                    all_beta(step, 0) = carDataset.dt * step;
                    last_beta         = val;
                    break;
                case 'r':
                    all_yawrate(step, 1) = val;
                    all_yawrate(step, 0) = carDataset.dt * step;
                    last_yawrate         = val;
                    break;
            };
        }
    };

    CarParameters cp;

    for (unsigned int timeStep = 0; timeStep < N; timeStep++, t += dt)
    {
        // sensor readings:
        const double u        = carDataset.vx[timeStep];
        const double delta    = carDataset.delta[timeStep];
        const double obs_r_k  = carDataset.yawRate[timeStep];
        const double obs_ay_k = carDataset.ay[timeStep];

        all_obs_r(timeStep, 0) = t;
        all_obs_r(timeStep, 1) = obs_r_k;

        const bool isLast = (timeStep == (N - 1));

        if (!isLast)
        {
            // ======= FactorBeta ======
            wholeFG.emplace_shared<FactorBeta>(
                cp, carDataset.dt, u, delta, noise_factor_beta, sBeta(timeStep),
                sBeta(timeStep + 1), sYawRate(timeStep));

            // ======= FactorRate ======
            wholeFG.emplace_shared<FactorRate>(
                cp, carDataset.dt, u, delta, noise_factor_rate,
                sYawRate(timeStep), sYawRate(timeStep + 1), sBeta(timeStep));

            wholeValues.insert(sBeta(timeStep + 1), last_beta);
            wholeValues.insert(sYawRate(timeStep + 1), last_yawrate);
        }

        // ======= YawRate observation======
        wholeFG.emplace_shared<gtsam::PriorFactor<double>>(
            sYawRate(timeStep), obs_r_k, noise_observation_rate);

        // ======= FactorLateralAcceleration ======
        wholeFG.emplace_shared<FactorLateralAcceleration>(
            cp, u, delta, obs_ay_k, noise_factor_ay, sBeta(timeStep),
            sYawRate(timeStep));

    }  // end for each time step

    tle1.stop();

    // ========= OPTIMIZE ==========
    gtsam::Values estimated;

    auto optParam          = gtsam::GaussNewtonParams();
    optParam.maxIterations = 1;

    const int fixedLagWindowLen = arg_fixed_lag_window.getValue();
    if (fixedLagWindowLen <= 0)
    {
        mrpt::system::CTimeLoggerEntry tle2(profiler, "2.optimize_fg");

        optParam.iterationHook =
            [&wholeFG](size_t iter, double errBef, double errAfter) {
                const auto N = wholeFG.size();
                std::cout << "LM iter #" << iter
                          << " rmse: " << std::sqrt(errBef / N) << " -> "
                          << std::sqrt(errAfter / N) << std::endl;
            };

        gtsam::GaussNewtonOptimizer optimizer(wholeFG, wholeValues, optParam);

        estimated = optimizer.optimize();
        lambda_Values_toAll(estimated);
    }
    else
    {
        mrpt::system::CTimeLoggerEntry tle2(
            profiler, "2.fixed_lag_optimize_fg");

        estimated = wholeValues;

        // Sort factors by time index:
        std::map<int, std::vector<gtsam::NonlinearFactor::shared_ptr>>
            time2factors;

        for (const auto& f : wholeFG)
        {
            std::set<int> ts;
            for (const auto& k : f->keys())
            {
                const auto idx = gtsam::Symbol(k).index();
                ts.insert(idx);
            }

            // all keys within window?
            if (*ts.rbegin() > (*ts.begin() + fixedLagWindowLen))
                continue;  // skip:

            time2factors[*ts.begin()].push_back(f);
        }

        for (int i = 0; (i + fixedLagWindowLen) < int(N); i++)
        {
            gtsam::NonlinearFactorGraph winFG;
            for (int off = 0; off < fixedLagWindowLen; off++)
                for (const auto& f : time2factors.at(i + off)) winFG += f;

            gtsam::Values winValues;

            for (int idx = i; idx <= (i + fixedLagWindowLen); idx++)
            {
                winValues.insert(sBeta(idx), wholeValues.at(sBeta(idx)));
                winValues.insert(sYawRate(idx), wholeValues.at(sYawRate(idx)));
            }

            // Add strong prior for new "initial" values:
            winFG.emplace_shared<gtsam::PriorFactor<double>>(
                sBeta(i), estimated.at<double>(sBeta(i)), noise_prior_window);

            winFG.emplace_shared<gtsam::PriorFactor<double>>(
                sYawRate(i), estimated.at<double>(sYawRate(i)),
                noise_prior_window);

            // std::cout << "\n ======== i=" << i << "\n";

#if 0
            optParam.iterationHook =
                [&winFG](size_t iter, double errBef, double errAfter) {
                    const auto N = winFG.size();
                    std::cout << "LM iter #" << iter
                              << " rmse: " << std::sqrt(errBef / N) << " -> "
                              << std::sqrt(errAfter / N) << std::endl;
                };
#endif
            gtsam::GaussNewtonOptimizer optimizer(winFG, winValues, optParam);

            const auto& winEstimated = optimizer.optimize();

#if 0
            winValues.print("==== INITIAL VALUES =====");
            winFG.print("==== WINDOW FG =====");
            winEstimated.print("==== FINAL VALUES =====");
#endif

            // Move to "wholeValues":
            for (const auto& kv : winEstimated)
                wholeValues.update(kv.key, kv.value);
        }
        lambda_Values_toAll(wholeValues);
    }

#if 0
    estimated.print("========= estimated values ======\n");

  const double errorAfter = smoother.getFactors().error(estimated);
  const auto numFactors = smoother.getFactors().size();

  std::cout << "ErrorAfter=" << errorAfter
            << " RMSE=" << std::sqrt(errorAfter / numFactors)
            << " numFactors=" << numFactors << "\n";
#endif

    auto lmbdPrintErr = [](const gtsam::Factor* /*factor*/,
                           double whitenedError, size_t /*index*/) -> bool {
        return true;  // whitenedError > 1e-3;
    };

    if (0)  // arg_show_factor_errors.isSet())
    {
        std::cout << "======== FACTOR ERRORS ============\n";
        // smoother.getFactors().
        wholeFG.printErrors(
            estimated, "FACTOR ERRORS: ", gtsam::DefaultKeyFormatter,
            lmbdPrintErr);
    }

    const auto prefix = arg_output_prefix.getValue();

    std::cout << "Saving results to TXT files with prefix '" << prefix
              << "'...\n";

    using namespace std::string_literals;  // "s"

    all_beta.saveToTextFile(
        prefix + "estimated_beta.txt"s, {}, false, "% TIMESTAMP beta\n");
    all_yawrate.saveToTextFile(
        prefix + "estimated_yawrate.txt"s, {}, false, "% TIMESTAMP r_estim\n");

    all_obs_r.saveToTextFile(
        prefix + "obs_yawrate.txt"s, {}, false, "% TIMESTAMP r_obs\n");

    // Evaluate beta estimation RMSE vs ground truth, if provided:
    if (!carDataset.beta_true.empty())
    {
        double err = 0;
        for (size_t i = 0; i < N; i++)
            err += mrpt::square(carDataset.beta_true[i] - all_beta(i, 1));
        if (N > 0) err /= N;
        const double beta_rmse = std::sqrt(err);
        printf(
            "Beta estimation RMSE vs ground truth=%16.6f rad (%16.6f deg)\n",
            beta_rmse, mrpt::RAD2DEG(beta_rmse));
    }
}

void test_factor_beta();
void test_factor_rate();
void test_factor_lateral_acceleration();

int main(int argc, char** argv)
{
    try
    {
        // -------------------------------------------
        // TODO: Move tests to independent files!
        test_factor_beta();
        test_factor_rate();
        test_factor_lateral_acceleration();
        // -------------------------------------------

        // Parse arguments:
        if (!cmd.parse(argc, argv))
            throw std::runtime_error("");  // should exit.

        slip_angle_smoother();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << mrpt::exception_to_str(e);
    }
}

void test_factor_beta()
{
// Use MRPT macros for testing as "EXPECT"
#define EXPECT ASSERT_
    const std::string name_ = "FactorBeta";

    const auto sBeta    = gtsam::SymbolGenerator('b');
    const auto sYawRate = gtsam::SymbolGenerator('r');

    auto noise = gtsam::noiseModel::Diagonal::Sigmas(gtsam::Vector1(1.0));
    CarParameters cp;
    const double  dt    = 1e-3;  // [s]
    const double  u     = 5.0;  // m/s
    const double  delta = mrpt::DEG2RAD(10.0);  // steering angle

    FactorBeta factor(cp, dt, u, delta, noise, sBeta(0), sBeta(1), sYawRate(0));

    // Set the linearization point
    gtsam::Values values;

    values.insert(sBeta(0), 1e-3);
    values.insert(sBeta(1), 2e-3);
    values.insert(sYawRate(0), 0.1);

    EXPECT_CORRECT_FACTOR_JACOBIANS(
        factor, values, 1e-7 /*diff*/, 1e-10 /*tolerance*/);
}

void test_factor_rate()
{
// Use MRPT macros for testing as "EXPECT"
#define EXPECT ASSERT_
    const std::string name_ = "FactorRate";

    const auto sBeta    = gtsam::SymbolGenerator('b');
    const auto sYawRate = gtsam::SymbolGenerator('r');

    auto noise = gtsam::noiseModel::Diagonal::Sigmas(gtsam::Vector1(1.0));
    CarParameters cp;
    const double  dt    = 1e-3;  // [s]
    const double  u     = 5.0;  // m/s
    const double  delta = mrpt::DEG2RAD(10.0);  // steering angle

    FactorRate factor(
        cp, dt, u, delta, noise, sYawRate(0), sYawRate(1), sBeta(0));

    // Set the linearization point
    gtsam::Values values;

    values.insert(sYawRate(0), 0.1);
    values.insert(sYawRate(1), 0.11);
    values.insert(sBeta(0), 0.01);

    EXPECT_CORRECT_FACTOR_JACOBIANS(
        factor, values, 1e-7 /*diff*/, 1e-10 /*tolerance*/);
}

void test_factor_lateral_acceleration()
{
// Use MRPT macros for testing as "EXPECT"
#define EXPECT ASSERT_
    const std::string name_ = "FactorLateralAcceleration";

    const auto sBeta    = gtsam::SymbolGenerator('b');
    const auto sYawRate = gtsam::SymbolGenerator('r');

    auto noise = gtsam::noiseModel::Diagonal::Sigmas(gtsam::Vector1(1.0));
    CarParameters cp;
    const double  ay    = 8.6;  // [m/s2]
    const double  u     = 5.0;  // m/s
    const double  delta = mrpt::DEG2RAD(10.0);  // steering angle

    FactorLateralAcceleration factor(
        cp, u, delta, ay, noise, sBeta(0), sYawRate(0));

    // Set the linearization point
    gtsam::Values values;

    values.insert(sBeta(0), 0.01);
    values.insert(sYawRate(0), 0.1);

    EXPECT_CORRECT_FACTOR_JACOBIANS(
        factor, values, 1e-7 /*diff*/, 1e-6 /*tolerance*/);
}
