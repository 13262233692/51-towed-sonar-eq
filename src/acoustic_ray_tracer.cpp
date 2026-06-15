#define _USE_MATH_DEFINES
#include "acoustic_ray_tracer.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace sonar {

AcousticRayTracer::AcousticRayTracer()
    : surface_depth_m_(0.0), bottom_depth_m_(5000.0),
      surface_reflectivity_(0.95), bottom_reflectivity_(0.5) {
    ssp_.depth_m = {0.0, 50.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0};
    ssp_.sound_speed_mps = {1530.0, 1510.0, 1490.0, 1480.0, 1485.0, 1500.0, 1520.0};
}

AcousticRayTracer::AcousticRayTracer(const SSProfile& ssp)
    : ssp_(ssp), surface_depth_m_(0.0), bottom_depth_m_(5000.0),
      surface_reflectivity_(0.95), bottom_reflectivity_(0.5) {}

void AcousticRayTracer::set_sound_speed_profile(const SSProfile& ssp) {
    ssp_ = ssp;
}

void AcousticRayTracer::set_environment(double surface_depth_m, double bottom_depth_m,
                                        double surface_reflectivity, double bottom_reflectivity) {
    surface_depth_m_ = surface_depth_m;
    bottom_depth_m_ = bottom_depth_m;
    surface_reflectivity_ = surface_reflectivity;
    bottom_reflectivity_ = bottom_reflectivity;
}

double AcousticRayTracer::interpolate_sound_speed(double depth_m) const {
    if (ssp_.depth_m.empty()) return SPEED_OF_SOUND_WATER;
    if (depth_m <= ssp_.depth_m.front()) return ssp_.sound_speed_mps.front();
    if (depth_m >= ssp_.depth_m.back()) return ssp_.sound_speed_mps.back();
    for (size_t i = 0; i + 1 < ssp_.depth_m.size(); ++i) {
        if (depth_m >= ssp_.depth_m[i] && depth_m <= ssp_.depth_m[i + 1]) {
            double t = (depth_m - ssp_.depth_m[i]) / (ssp_.depth_m[i + 1] - ssp_.depth_m[i]);
            return ssp_.sound_speed_mps[i] * (1 - t) + ssp_.sound_speed_mps[i + 1] * t;
        }
    }
    return SPEED_OF_SOUND_WATER;
}

double AcousticRayTracer::snell_refraction_angle(double theta1, double c1, double c2) {
    double sin_theta2 = (c2 / c1) * std::sin(theta1);
    if (sin_theta2 >= 1.0) return M_PI / 2.0;
    if (sin_theta2 <= -1.0) return -M_PI / 2.0;
    return std::asin(sin_theta2);
}

Vec3 AcousticRayTracer::refract_direction(const Vec3& dir, double c1, double c2,
                                           const Vec3& normal) {
    double cos_theta1 = -dir.dot(normal);
    cos_theta1 = std::max(-1.0, std::min(1.0, cos_theta1));
    double theta1 = std::acos(cos_theta1);
    double theta2 = snell_refraction_angle(theta1, c1, c2);
    Vec3 tangent = dir + normal * cos_theta1;
    double tlen = tangent.norm();
    if (tlen < 1e-12) return normal * (-1.0);
    tangent = tangent * (1.0 / tlen);
    return tangent * std::sin(theta2) - normal * std::cos(theta2);
}

AcousticRay AcousticRayTracer::trace_single_ray(const PlatformState& state,
                                                 const Vec3& launch_pt,
                                                 const Vec3& launch_dir_raw,
                                                 double max_bounces) {
    AcousticRay ray;
    ray.launch_point = launch_pt;
    double dnorm = launch_dir_raw.norm();
    Vec3 dir = (dnorm > 1e-12) ? launch_dir_raw * (1.0 / dnorm) : Vec3(0, 0, -1);
    ray.launch_direction = dir;

    Vec3 pos = launch_pt;
    double total_dist = 0.0;
    double total_time = 0.0;
    size_t bounces_s = 0, bounces_b = 0;
    double loss_db = 0.0;
    Vec3 normal_surface(0, 0, -1);
    Vec3 normal_bottom(0, 0, 1);
    const double step_m = 0.5;
    const size_t max_steps = 200000;
    const double target_hit_radius_m = 50.0;
    double grazing_angle = 0.0;
    double incident_angle = 0.0;

    for (size_t step = 0; step < max_steps; ++step) {
        double cur_c = interpolate_sound_speed(pos.z);
        Vec3 next_pos = pos + dir * step_m;
        double dist_inc = step_m;

        if (next_pos.z <= surface_depth_m_) {
            double t_hit = (surface_depth_m_ - pos.z) / (next_pos.z - pos.z);
            t_hit = std::max(0.0, std::min(1.0, t_hit));
            Vec3 hit_pt = pos + (next_pos - pos) * t_hit;
            dist_inc = (hit_pt - pos).norm();
            total_dist += dist_inc;
            total_time += dist_inc / cur_c;
            grazing_angle = std::acos(std::max(-1.0, std::min(1.0, -dir.z)));
            pos = hit_pt;
            dir = dir - normal_surface * 2.0 * dir.dot(normal_surface);
            loss_db += -20.0 * std::log10(std::max(1e-6, surface_reflectivity_));
            bounces_s++;
            if (static_cast<double>(bounces_s + bounces_b) > max_bounces) break;
            continue;
        }
        if (next_pos.z >= bottom_depth_m_) {
            double t_hit = (bottom_depth_m_ - pos.z) / (next_pos.z - pos.z);
            t_hit = std::max(0.0, std::min(1.0, t_hit));
            Vec3 hit_pt = pos + (next_pos - pos) * t_hit;
            dist_inc = (hit_pt - pos).norm();
            total_dist += dist_inc;
            total_time += dist_inc / cur_c;
            grazing_angle = std::acos(std::max(-1.0, std::min(1.0, dir.z)));
            pos = hit_pt;
            dir = dir - normal_bottom * 2.0 * dir.dot(normal_bottom);
            loss_db += -20.0 * std::log10(std::max(1e-6, bottom_reflectivity_));
            bounces_b++;
            if (static_cast<double>(bounces_s + bounces_b) > max_bounces) break;
            continue;
        }

        double next_c = interpolate_sound_speed(next_pos.z);
        Vec3 grad_dir = (dir.z >= 0) ? normal_bottom : normal_surface;
        dir = refract_direction(dir, cur_c, next_c, grad_dir);
        double dn = dir.norm();
        if (dn > 1e-12) dir = dir * (1.0 / dn);

        total_dist += dist_inc;
        total_time += dist_inc / cur_c;
        pos = next_pos;

        Vec3 to_target = state.target_position_m - pos;
        if (to_target.norm() < target_hit_radius_m) {
            Vec3 rel_vel = state.target_velocity_mps - state.sub_velocity_mps;
            Vec3 rdir = pos - state.sub_position_m;
            double rd = rdir.norm();
            if (rd > 1e-6) rdir = rdir * (1.0 / rd);
            ray.radial_velocity_mps = rel_vel.dot(rdir);
            ray.doppler_shift_hz = ray.radial_velocity_mps / SPEED_OF_SOUND_WATER * state.carrier_frequency_hz;
            incident_angle = std::acos(std::max(-1.0, std::min(1.0,
                -dir.dot(rdir))));
            break;
        }
    }

    if (total_dist > 0 && ray.radial_velocity_mps == 0.0 && ray.doppler_shift_hz == 0.0) {
        Vec3 rel_vel = state.target_velocity_mps - state.sub_velocity_mps;
        Vec3 rdir = pos - state.sub_position_m;
        double rd = rdir.norm();
        if (rd > 1e-6) {
            rdir = rdir * (1.0 / rd);
            ray.radial_velocity_mps = rel_vel.dot(rdir);
            ray.doppler_shift_hz = ray.radial_velocity_mps / SPEED_OF_SOUND_WATER * state.carrier_frequency_hz;
        }
    }

    ray.travel_time_s = total_time;
    ray.path_length_m = total_dist;
    ray.grazing_angle_rad = grazing_angle;
    ray.incident_angle_rad = incident_angle;
    ray.reflection_loss_db = loss_db;
    ray.surface_bounces = bounces_s;
    ray.bottom_bounces = bounces_b;
    double amp = std::pow(10.0, -loss_db / 20.0) / std::max(1.0, total_dist);
    ray.complex_weight = amp;
    return ray;
}

std::vector<AcousticRay> AcousticRayTracer::trace_rays(const PlatformState& state,
                                                        const ArrayDeformation& array,
                                                        size_t num_rays_per_quadrant,
                                                        double max_bounces) {
    std::vector<AcousticRay> rays;
    Vec3 launch_base = state.sub_position_m;
    if (!array.element_positions_m.empty()) {
        launch_base = array.element_positions_m[array.element_positions_m.size() / 2];
    }
    Vec3 to_tgt = state.target_position_m - launch_base;
    double tgt_dist = to_tgt.norm();
    if (tgt_dist < 1e-6) return rays;
    Vec3 forward = to_tgt * (1.0 / tgt_dist);
    Vec3 up(0, 0, -1);
    Vec3 right = Vec3(forward.y * up.z - forward.z * up.y,
                       forward.z * up.x - forward.x * up.z,
                       forward.x * up.y - forward.y * up.x);
    double rl = right.norm();
    if (rl < 1e-6) { right = Vec3(1, 0, 0); } else { right = right * (1.0 / rl); }
    Vec3 real_up = Vec3(right.y * forward.z - right.z * forward.y,
                         right.z * forward.x - right.x * forward.z,
                         right.x * forward.y - right.y * forward.x);

    size_t id = 0;
    double dtheta = M_PI_2 / std::max((size_t)1, num_rays_per_quadrant);
    for (int az = -(int)num_rays_per_quadrant; az <= (int)num_rays_per_quadrant; ++az) {
        double azi = az * dtheta * 0.5;
        for (int el = -(int)num_rays_per_quadrant; el <= (int)num_rays_per_quadrant; ++el) {
            double ele = el * dtheta * 0.25;
            Vec3 dir = forward * std::cos(ele) * std::cos(azi)
                      + right * std::cos(ele) * std::sin(azi)
                      + real_up * std::sin(ele);
            AcousticRay r = trace_single_ray(state, launch_base, dir, max_bounces);
            r.ray_id = id++;
            rays.push_back(r);
        }
    }
    return rays;
}

std::vector<double> AcousticRayTracer::compute_radial_velocities(
    const PlatformState& state, const std::vector<AcousticRay>& rays) {
    std::vector<double> result;
    result.reserve(rays.size());
    for (const auto& r : rays) result.push_back(r.radial_velocity_mps);
    return result;
}

std::vector<double> AcousticRayTracer::compute_doppler_shifts(
    const PlatformState& state, const std::vector<AcousticRay>& rays) {
    std::vector<double> result;
    result.reserve(rays.size());
    for (const auto& r : rays) result.push_back(r.doppler_shift_hz);
    return result;
}

}
