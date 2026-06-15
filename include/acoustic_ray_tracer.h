#ifndef ACOUSTIC_RAY_TRACER_H
#define ACOUSTIC_RAY_TRACER_H

#include "sonar_types.h"
#include <vector>

namespace sonar {

class AcousticRayTracer {
public:
    AcousticRayTracer();
    explicit AcousticRayTracer(const SSProfile& ssp);

    void set_sound_speed_profile(const SSProfile& ssp);
    void set_environment(double surface_depth_m, double bottom_depth_m,
                         double surface_reflectivity, double bottom_reflectivity);

    std::vector<AcousticRay> trace_rays(const PlatformState& state,
                                        const ArrayDeformation& array,
                                        size_t num_rays_per_quadrant = 16,
                                        double max_bounces = 2);

    std::vector<double> compute_radial_velocities(const PlatformState& state,
                                                   const std::vector<AcousticRay>& rays);

    std::vector<double> compute_doppler_shifts(const PlatformState& state,
                                                const std::vector<AcousticRay>& rays);

    double interpolate_sound_speed(double depth_m) const;

private:
    SSProfile ssp_;
    double surface_depth_m_;
    double bottom_depth_m_;
    double surface_reflectivity_;
    double bottom_reflectivity_;

    AcousticRay trace_single_ray(const PlatformState& state,
                                  const Vec3& launch_pt,
                                  const Vec3& launch_dir,
                                  double max_bounces);

    Vec3 refract_direction(const Vec3& dir, double c1, double c2,
                           const Vec3& normal);

    double snell_refraction_angle(double theta1, double c1, double c2);
};

}

#endif
