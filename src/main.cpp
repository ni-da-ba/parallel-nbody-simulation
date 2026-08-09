#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "distribution.hpp"
#include "octree.hpp"
#include "particle.hpp"
#include "particle_io.hpp"
#include "simulation.hpp"
#include "vec3.hpp"

namespace {
  constexpr double G = 1.0;
  constexpr double SOFTENING = 1e-9;

  std::vector<Particle> make_initial_particles(int n) {
    std::vector<Particle> particles;
    particles.reserve(n);

    if (n == 3) {
      particles.emplace_back(Vec3(-1.0, 0.0, 0.0), Vec3(0.0, 0.2, 0.0), 1.0);
      particles.emplace_back(Vec3( 1.0, 0.0, 0.0), Vec3(0.0,-0.2, 0.0), 1.0);
      particles.emplace_back(Vec3( 0.0, 1.0, 0.0), Vec3(-0.2, 0.0, 0.0), 0.5);
      return particles;
    }

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> pos_dist(-1.0, 1.0);
    std::uniform_real_distribution<double> vel_dist(-0.05, 0.05);
    std::uniform_real_distribution<double> mass_dist(0.5, 1.5);

    for (int i = 0; i < n; ++i) {
      Vec3 position(pos_dist(rng), pos_dist(rng), pos_dist(rng));
      Vec3 velocity(vel_dist(rng), vel_dist(rng), vel_dist(rng));
      double mass = mass_dist(rng);
      particles.emplace_back(position, velocity, mass);
    }

    return particles;
  }

  void print_usage_if_needed(int rank) {
    if (rank == 0) {
      std::cout
	<< "Usage: ./nbody [mode] [N] [steps] [dt] [theta] [snapshot_interval] [snapshot_prefix]\n"
	<< "  mode: allpairs | barneshut\n"
	<< "Defaults: mode=allpairs, N=3, steps=5, dt=0.01, theta=0.5, snapshot_interval=0, snapshot_prefix=snapshot\n\n";
    }
  }

  void save_snapshot_csv(
			 const std::vector<Particle>& particles,
			 int step,
			 const std::string& prefix
			 ) {
    std::ostringstream filename;
    filename << prefix << "_step" << std::setw(4) << std::setfill('0') << step << ".csv";

    std::ofstream out(filename.str());
    out << "id,x,y,z,vx,vy,vz,mass\n";

    for (std::size_t i = 0; i < particles.size(); ++i) {
      const auto& p = particles[i];
      out << i << ','
	  << p.position.x << ','
	  << p.position.y << ','
	  << p.position.z << ','
	  << p.velocity.x << ','
	  << p.velocity.y << ','
	  << p.velocity.z << ','
	  << p.mass << '\n';
    }
  }

  void step_particles_barnes_hut_from_global(
					     std::vector<Particle>& local_particles,
					     const std::vector<Particle>& global_particles,
					     int global_start_index,
					     double dt,
        double theta
					     ) {
    Cube root_region = make_bounding_cube(global_particles);
    OctreeNode root(root_region);

    for (int i = 0; i < static_cast<int>(global_particles.size()); ++i) {
      root.insert(i, global_particles);
    }

    root.compute_mass_properties(global_particles);

    std::vector<Vec3> accelerations(local_particles.size(), Vec3(0.0, 0.0, 0.0));

    for (std::size_t i = 0; i < local_particles.size(); ++i) {
      const int global_index = global_start_index + static_cast<int>(i);

      accelerations[i] = root.compute_acceleration_on_particle(
							       global_index, global_particles, theta, SOFTENING, G
							       );
    }

    for (std::size_t i = 0; i < local_particles.size(); ++i) {
      local_particles[i].acceleration = accelerations[i];
      local_particles[i].velocity += local_particles[i].acceleration * dt;
      local_particles[i].position += local_particles[i].velocity * dt;
    }
  }

  std::vector<Particle> gather_global_particles(
						const std::vector<Particle>& local_particles,
						int global_particle_count,
						const std::vector<int>& packed_counts,
						const std::vector<int>& packed_offsets
						) {
    std::vector<double> packed_local = pack_particles(local_particles);
    std::vector<double> packed_global(global_particle_count * PARTICLE_DOUBLES);

    MPI_Allgatherv(
		   packed_local.data(),
		   static_cast<int>(packed_local.size()),
		   MPI_DOUBLE,
		   packed_global.data(),
		   packed_counts.data(),
		   packed_offsets.data(),
		   MPI_DOUBLE,
            MPI_COMM_WORLD
		   );

    return unpack_particles(packed_global);
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  int size = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::string mode = "allpairs";
  int global_particle_count = 3;
  int steps = 5;
  double dt = 0.01;
  double theta = 0.5;
  int snapshot_interval = 0;
  std::string snapshot_prefix = "snapshot";

  if (argc >= 2) {
    mode = argv[1];
  }
  if (argc >= 3) {
    global_particle_count = std::atoi(argv[2]);
  }
  if (argc >= 4) {
    steps = std::atoi(argv[3]);
  }
  if (argc >= 5) {
    dt = std::atof(argv[4]);
  }
  if (argc >= 6) {
    theta = std::atof(argv[5]);
  }
  if (argc >= 7) {
    snapshot_interval = std::atoi(argv[6]);
  }
  if (argc >= 8) {
    snapshot_prefix = argv[7];
  }

  if ((mode != "allpairs" && mode != "barneshut") ||
        global_particle_count <= 0 || steps < 0 || dt <= 0.0 ||
      theta <= 0.0 || snapshot_interval < 0) {
    print_usage_if_needed(rank);
    if (rank == 0) {
      std::cerr << "Error: mode must be allpairs or barneshut; require N > 0, steps >= 0, dt > 0, theta > 0, snapshot_interval >= 0.\n";
    }
    MPI_Finalize();
    return 1;
  }

  std::vector<Particle> global_particles;
  std::vector<double> packed_global_particles;

  if (rank == 0) {
    global_particles = make_initial_particles(global_particle_count);
    packed_global_particles = pack_particles(global_particles);
  }

  MPI_Bcast(&global_particle_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> particle_counts = all_local_counts(global_particle_count, size);
  std::vector<int> particle_offsets = all_local_offsets(global_particle_count, size);

  std::vector<int> packed_counts(size);
  std::vector<int> packed_offsets(size);

  for (int r = 0; r < size; ++r) {
    packed_counts[r] = particle_counts[r] * PARTICLE_DOUBLES;
    packed_offsets[r] = particle_offsets[r] * PARTICLE_DOUBLES;
  }

  const int local_particle_n = particle_counts[rank];
  std::vector<double> packed_local_particles(local_particle_n * PARTICLE_DOUBLES);

  MPI_Scatterv(
	       rank == 0 ? packed_global_particles.data() : nullptr,
	       packed_counts.data(),
	       packed_offsets.data(),
	       MPI_DOUBLE,
	       packed_local_particles.data(),
	       static_cast<int>(packed_local_particles.size()),
	       MPI_DOUBLE,
	       0,
        MPI_COMM_WORLD
	       );

  std::vector<Particle> local_particles = unpack_particles(packed_local_particles);

  if (rank == 0) {
    std::cout << "MPI initialized with " << size << " ranks.\n";
    std::cout << "Running mode = " << mode
	      << ", N = " << global_particle_count
	      << ", steps = " << steps
	      << ", dt = " << dt
	      << ", theta = " << theta
	      << ", snapshot_interval = " << snapshot_interval
	      << ", snapshot_prefix = " << snapshot_prefix
	      << "\n\n";
  }

  Vec3 initial_momentum(0.0, 0.0, 0.0);
  double initial_mass = 0.0;

  global_particles = gather_global_particles(
					     local_particles, global_particle_count, packed_counts, packed_offsets
					     );

  if (rank == 0) {
    initial_mass = compute_total_mass(global_particles);
    initial_momentum = compute_total_momentum(global_particles);

    if (global_particle_count <= 10) {
      print_particles(global_particles, 0);
    } else {
      std::cout << "Initial ";
      print_system_statistics(global_particles, 0);
      std::cout << '\n';
    }

    if (snapshot_interval > 0) {
      save_snapshot_csv(global_particles, 0, snapshot_prefix);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  const double start_time = MPI_Wtime();

  for (int step = 1; step <= steps; ++step) {
    global_particles = gather_global_particles(
					       local_particles, global_particle_count, packed_counts, packed_offsets
					       );

    if (mode == "allpairs") {
      step_particles_from_global(local_particles, global_particles, dt);
    } else {
      step_particles_barnes_hut_from_global(
					    local_particles,
					    global_particles,
					    particle_offsets[rank],
					    dt,
                theta
					    );
    }

    const bool need_snapshot = (snapshot_interval > 0 && step % snapshot_interval == 0);
    const bool need_small_print = (global_particle_count <= 10);

    if (need_snapshot || need_small_print) {
      std::vector<Particle> updated_global_particles = gather_global_particles(
									       local_particles, global_particle_count, packed_counts, packed_offsets
									       );

      if (rank == 0) {
	if (need_small_print) {
	  print_particles(updated_global_particles, step);
	}
	if (need_snapshot) {
	  save_snapshot_csv(updated_global_particles, step, snapshot_prefix);
	}
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  const double end_time = MPI_Wtime();
  const double local_elapsed = end_time - start_time;

  double max_elapsed = 0.0;
  MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  global_particles = gather_global_particles(
					     local_particles, global_particle_count, packed_counts, packed_offsets
					     );

  if (rank == 0) {
    const double final_mass = compute_total_mass(global_particles);
    const Vec3 final_momentum = compute_total_momentum(global_particles);

    const double mass_change = final_mass - initial_mass;
    const Vec3 momentum_change(
			       final_momentum.x - initial_momentum.x,
			       final_momentum.y - initial_momentum.y,
            final_momentum.z - initial_momentum.z
			       );

    const double p_mag_change = std::sqrt(
            momentum_change.x * momentum_change.x +
            momentum_change.y * momentum_change.y +
            momentum_change.z * momentum_change.z
					  );

    if (global_particle_count > 10) {
      std::cout << "Final ";
      print_system_statistics(global_particles, steps);
      std::cout << '\n';
    }

    std::cout << "Mass change: " << mass_change << '\n';
    std::cout << "Momentum change: ("
	      << momentum_change.x << ", "
	      << momentum_change.y << ", "
	      << momentum_change.z << ")\n";
    std::cout << "Momentum change magnitude: " << p_mag_change << "\n\n";

    std::cout << "Elapsed time: " << max_elapsed << " seconds\n";
    if (steps > 0) {
      std::cout << "Average time per step: " << (max_elapsed / steps) << " seconds\n";
    }

    std::cout
      << "SUMMARY "
      << "mode=" << mode
      << " ranks=" << size
      << " N=" << global_particle_count
      << " steps=" << steps
      << " dt=" << dt
      << " theta=" << theta
      << " elapsed=" << max_elapsed
      << " avg_step=" << (steps > 0 ? max_elapsed / steps : 0.0)
      << " mass_change=" << mass_change
      << " px_change=" << momentum_change.x
      << " py_change=" << momentum_change.y
      << " pz_change=" << momentum_change.z
      << " p_mag_change=" << p_mag_change
      << '\n';
  }

  MPI_Finalize();
  return 0;
}
