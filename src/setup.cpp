#include "setup.hpp"
#include <fstream>
#include <vector>
#include <cmath>
using std::vector;

#ifdef VREMAN_SGS
// Helper function to compute Vreman SGS viscosity at a point from velocity field
// Returns nu_sgs in LBM units
inline float compute_vreman_nu_sgs(const LBM& lbm, uint x, uint y, uint z) {
	const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();

	// Boundary check - need neighbors
	if(x < 1 || x >= Nx-1 || y < 1 || y >= Ny-1 || z < 1 || z >= Nz-1) return 0.0f;

	// Neighbor indices
	auto idx = [&](uint xx, uint yy, uint zz) -> ulong {
		return (ulong)xx + (ulong)yy * (ulong)Nx + (ulong)zz * (ulong)Nx * (ulong)Ny;
	};

	// Read neighbor velocities
	const float ux_xp = lbm.u.x[idx(x+1,y,z)], ux_xm = lbm.u.x[idx(x-1,y,z)];
	const float ux_yp = lbm.u.x[idx(x,y+1,z)], ux_ym = lbm.u.x[idx(x,y-1,z)];
	const float ux_zp = lbm.u.x[idx(x,y,z+1)], ux_zm = lbm.u.x[idx(x,y,z-1)];

	const float uy_xp = lbm.u.y[idx(x+1,y,z)], uy_xm = lbm.u.y[idx(x-1,y,z)];
	const float uy_yp = lbm.u.y[idx(x,y+1,z)], uy_ym = lbm.u.y[idx(x,y-1,z)];
	const float uy_zp = lbm.u.y[idx(x,y,z+1)], uy_zm = lbm.u.y[idx(x,y,z-1)];

	const float uz_xp = lbm.u.z[idx(x+1,y,z)], uz_xm = lbm.u.z[idx(x-1,y,z)];
	const float uz_yp = lbm.u.z[idx(x,y+1,z)], uz_ym = lbm.u.z[idx(x,y-1,z)];
	const float uz_zp = lbm.u.z[idx(x,y,z+1)], uz_zm = lbm.u.z[idx(x,y,z-1)];

	// Velocity gradients (central differences)
	const float a11 = 0.5f*(ux_xp - ux_xm); // du/dx
	const float a12 = 0.5f*(ux_yp - ux_ym); // du/dy
	const float a13 = 0.5f*(ux_zp - ux_zm); // du/dz
	const float a21 = 0.5f*(uy_xp - uy_xm); // dv/dx
	const float a22 = 0.5f*(uy_yp - uy_ym); // dv/dy
	const float a23 = 0.5f*(uy_zp - uy_zm); // dv/dz
	const float a31 = 0.5f*(uz_xp - uz_xm); // dw/dx
	const float a32 = 0.5f*(uz_yp - uz_ym); // dw/dy
	const float a33 = 0.5f*(uz_zp - uz_zm); // dw/dz

	// Alpha invariant
	const float alpha_sq = a11*a11 + a12*a12 + a13*a13
	                     + a21*a21 + a22*a22 + a23*a23
	                     + a31*a31 + a32*a32 + a33*a33;

	if(alpha_sq < 1e-12f) return 0.0f;

	// Beta tensor
	const float b11 = a11*a11 + a21*a21 + a31*a31;
	const float b22 = a12*a12 + a22*a22 + a32*a32;
	const float b33 = a13*a13 + a23*a23 + a33*a33;
	const float b12 = a11*a12 + a21*a22 + a31*a32;
	const float b13 = a11*a13 + a21*a23 + a31*a33;
	const float b23 = a12*a13 + a22*a23 + a32*a33;

	// B_beta invariant
	const float B_beta = b11*b22 - b12*b12 + b11*b33 - b13*b13 + b22*b33 - b23*b23;

	float nu_sgs = 0.0f;
	if(B_beta > 1e-12f) {
		// Vreman constant (same as in kernel.cpp def_Cv)
		const float Cv = 0.07f;
		nu_sgs = Cv * std::sqrt(B_beta / alpha_sq);
	}

	// Stability floor (same as in kernel.cpp def_Cv_floor)
	const float Cv_floor = 1e-5f;
	return std::max(nu_sgs, Cv_floor);
}
#endif // VREMAN_SGS


#ifdef BENCHMARK
#include "info.hpp"
void main_setup() { // benchmark; required extensions in defines.hpp: BENCHMARK, optionally FP16S or FP16C
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	uint mlups = 0u; {

		//LBM lbm( 32u,  32u,  32u, 1.0f);
		//LBM lbm( 64u,  64u,  64u, 1.0f);
		//LBM lbm(128u, 128u, 128u, 1.0f);
		LBM lbm(256u, 256u, 256u, 1.0f); // default
		//LBM lbm(384u, 384u, 384u, 1.0f);
		//LBM lbm(512u, 512u, 512u, 1.0f);

		//const uint memory = 1488u; // memory occupation in MB (for multi-GPU benchmarks: make this close to as large as the GPU's VRAM capacity)
		//const uint3 lbm_N = (resolution(float3(1.0f, 1.0f, 1.0f), memory)/4u)*4u; // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
		//LBM lbm(1u*lbm_N.x, 1u*lbm_N.y, 1u*lbm_N.z, 1u, 1u, 1u, 1.0f); // 1 GPU
		//LBM lbm(2u*lbm_N.x, 1u*lbm_N.y, 1u*lbm_N.z, 2u, 1u, 1u, 1.0f); // 2 GPUs
		//LBM lbm(2u*lbm_N.x, 2u*lbm_N.y, 1u*lbm_N.z, 2u, 2u, 1u, 1.0f); // 4 GPUs
		//LBM lbm(2u*lbm_N.x, 2u*lbm_N.y, 2u*lbm_N.z, 2u, 2u, 2u, 1.0f); // 8 GPUs

		// #########################################################################################################################################################################################
		for(uint i=0u; i<1000u; i++) {
			lbm.run(10u, 1000u*10u);
			mlups = max(mlups, to_uint((double)lbm.get_N()*1E-6/info.runtime_lbm_timestep_smooth));
		}
	} // make lbm object go out of scope to free its memory
	print_info("Peak MLUPs/s = "+to_string(mlups));
#if defined(_WIN32)
	wait();
#endif // Windows
} /**/
#endif // BENCHMARK



#ifndef BENCHMARK
void main_setup() { // Simple rectangular jet test; required extensions: EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint memory = 6000u; // MB VRAM (increased for better resolution)
	// Plane jet with tanh inlet profile - no nozzle, no correction factor needed
	// U(y) = U_j/2 * (1 + tanh((h/2 - |y|) / (2*theta))), h/theta = 20
	const float lbm_u = 0.05f; // inlet velocity in LBM units (jet centerline velocity)

	// Domain: ~28h x 20h x 4h (spanwise-periodic thin slab)
	// Z reduced to 0.2 ratio for Lz/h ≈ 4 (periodic, per Stanley & Sarkar 2002)
	// Resolution gain: h_cells ~21 → ~37 (76% improvement)
	// Max usable x/h ≈ 25 (sufficient for transition + self-similar regions)
	const uint3 lbm_N = resolution(float3(1.4f, 1.0f, 0.2f), memory);

	// Re = U * h / nu => nu = U * h / Re
	// Compute viscosity using actual grid resolution (h_cells = Ny/20)
	// Using SUBGRID turbulence model for stability at high Re
	const float Re_target = 20100.0f;
	const uint h_cells_for_nu = lbm_N.y / 20u; // actual h_cells from resolved grid
	const float lbm_nu = lbm_u * (float)h_cells_for_nu / Re_target;

	LBM lbm(lbm_N.x, lbm_N.y, lbm_N.z, lbm_nu);

	// Plane jet dimensions (spanwise-periodic, no physical nozzle)
	const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	const uint h_cells = Ny / 20u;  // nozzle height = 1/20 of domain height

	const float Re_actual = lbm_u * (float)h_cells / lbm_nu;

	print_info("Grid: " + to_string(Nx) + " x " + to_string(Ny) + " x " + to_string(Nz));
	print_info("Plane jet: h = " + to_string(h_cells) + " cells, Lz/h = " + to_string(Nz / h_cells) + " (periodic)");
	print_info("Inlet: tanh profile, h/theta = 20, no physical nozzle");
	print_info("Inlet velocity = " + to_string(lbm_u));
	print_info("Viscosity = " + to_string(lbm_nu));
	print_info("Reynolds number = " + to_string((int)Re_actual));

	// ###################################################################################### define geometry ######################################################################################
	parallel_for(lbm.get_N(), [&](ulong n) {
		uint x = 0u, y = 0u, z = 0u;
		lbm.coordinates(n, x, y, z);

		// Center of domain (y only - z is periodic, no center needed)
		const int cy = (int)Ny / 2;

		// ============ Inlet BC (x = 0) - Tanh velocity profile with nozzle lip ============
		// U(y) = U_j/2 * (1 + tanh((h/2 - |y-cy|) / (2*theta))), h/theta = 20
		// Thin nozzle lip (TYPE_S) at y = cy ± h/2 to create sharp edge for K-H instability
		if(x == 0u) {
			const uint dy_abs = (uint)abs((int)y - cy);
			if(dy_abs == h_cells / 2u || dy_abs == h_cells / 2u + 1u) {
				// Nozzle lip walls (thin splitter plate at jet edges)
				lbm.flags[n] = TYPE_S;
			} else if(dy_abs < h_cells / 2u) {
				// Inside jet: tanh profile (core region)
				const float y_dist = fabs((float)y - (float)cy);
				const float theta = (float)h_cells / 20.0f;
				const float tanh_val = tanhf(((float)h_cells / 2.0f - y_dist) / (2.0f * theta));
				const float u_profile = lbm_u * 0.5f * (1.0f + tanh_val);
				lbm.flags[n] = TYPE_X; // DFM turbulent inlet
				lbm.u.x[n] = u_profile;
				lbm.u.y[n] = 0.0f;
				lbm.u.z[n] = 0.0f;
			} else {
				// Outside jet: ambient (quiescent)
				lbm.flags[n] = TYPE_E;
				lbm.u.x[n] = 0.0f;
				lbm.u.y[n] = 0.0f;
				lbm.u.z[n] = 0.0f;
			}
		}
		// ============ Outlet BC (last 3 cells) ============
		// FIX: Check outlet FIRST and unconditionally to ensure ALL outlet cells use TYPE_Y
		// This prevents lateral TYPE_E boundaries from creating a "zero-velocity frame" at outlet edges
		else if(x >= Nx - 3u) {
			lbm.flags[n] = TYPE_Y; // Convective outlet (Orlanski BC) - ALL cells including edges!
			// FIX: Initialize with ZERO velocity - let flow develop naturally from upstream
			// Setting a non-zero velocity creates artificial suction!
			lbm.u.x[n] = 0.0f;
			lbm.u.y[n] = 0.0f;
			lbm.u.z[n] = 0.0f;
		}
		// ============ Lateral boundaries - y only (z is periodic) ============
		else if(y <= 1u || y >= Ny-2u) {
			lbm.flags[n] = TYPE_E; // Ambient
			lbm.u.x[n] = 0.0f;
			lbm.u.y[n] = 0.0f;
			lbm.u.z[n] = 0.0f;
		}
		// ============ Interior fluid ============
		else {
			// Initialize with zero velocity
			lbm.u.x[n] = 0.0f;
			lbm.u.y[n] = 0.0f;
			lbm.u.z[n] = 0.0f;
		}
	});
	// ####################################################################### run simulation ##########################################################################
	lbm.graphics.visualization_modes = VIS_FIELD;

	// Calculate flow-through time (domain length / inlet velocity)
	const float domain_length = (float)Nx;
	const float flow_through_time = domain_length / lbm_u; // timesteps for one flow-through
	const uint sample_interval = 500u; // sample every 500 timesteps (>2x integral time scale for independent samples)
	const uint report_interval = (uint)(flow_through_time / 2.0f); // report every 0.5 flow-through

	// Timing parameters
	const float warmup_FT = 5.0f; // warmup: 5 flow-through times (ignore initial transient)
	const float averaging_FT = 30.0f; // averaging: 30 flow-through times (total 35 FT, then auto-stop)
	const ulong warmup_steps = (ulong)(warmup_FT * flow_through_time);
	const ulong total_steps = (ulong)((warmup_FT + averaging_FT) * flow_through_time);

	// Measurement locations
	const uint cy = Ny / 2u; // centerline y
	const uint cz = Nz / 2u; // centerline z

	// Spanwise averaging: average over ALL z-planes (z is periodic)
	// With periodic z BCs, all z-planes are statistically equivalent → maximum convergence
	const int z_avg_count = (int)Nz; // all z-planes

	// x/h stations for measurements (relative to inlet at x=0)
	// Paper Fig 4: x/h = 0, 2, 4, 6, 8, 10
	// Paper Fig 7: continuous x/h from 0 to ~30 with 0.2 increment
	const float xh_increment = 0.2f;
	const float xh_start = 0.0f;  // Start at inlet (no nozzle)
	const float xh_end = 30.0f;
	const int num_x_stations = (int)((xh_end - xh_start) / xh_increment) + 1; // x/h = 0, 0.2, ... , 30.0
	vector<uint> x_stations(num_x_stations);
	vector<float> xh_values(num_x_stations);
	for(int i = 0; i < num_x_stations; i++) {
		xh_values[i] = xh_start + (float)i * xh_increment;
		// x position in cells: x/h=0 is at x=0, so x = xh * h_cells
		const int x_pos = (int)round(xh_values[i] * (float)h_cells);
		x_stations[i] = (x_pos >= 0) ? (uint)x_pos : 0u;  // clamp to valid range
	}

	// Lateral profile extent (for Fig 4)
	const uint lateral_extent = 5u * h_cells; // measure from centerline to ±5h

	// ============ Statistics accumulators ============
	// For centerline (Fig 7a, 7b): store Uc at each x/h
	// With spanwise averaging, each timestep contributes z_avg_count samples
	vector<double> sum_Uc(num_x_stations, 0.0);
	vector<double> sum_Uc2(num_x_stations, 0.0);

	// For lateral profiles (Fig 4): store U at each (x, y), averaged over z-planes
	// We'll measure from y = cy - lateral_extent to y = cy + lateral_extent
	const uint lateral_points = 2u * lateral_extent + 1u;
	vector<vector<double>> sum_U_lateral(num_x_stations, vector<double>(lateral_points, 0.0));
	vector<vector<double>> sum_U2_lateral(num_x_stations, vector<double>(lateral_points, 0.0));

	// For volume flow (Fig 10): integrate U over cross-section at each x/h
	// Q is computed per z-plane, then averaged
	vector<double> sum_Q(num_x_stations, 0.0);

	ulong num_time_samples = 0ul; // number of timesteps sampled
	ulong num_total_samples = 0ul; // total samples = num_time_samples × z_avg_count

	print_info("Flow-through time = " + to_string((uint)flow_through_time) + " timesteps");
	print_info("Warmup: " + to_string(warmup_FT, 1u) + " FT (" + to_string(warmup_steps) + " steps)");
	print_info("Averaging: " + to_string(averaging_FT, 1u) + " FT");
	print_info("Total: " + to_string(warmup_FT + averaging_FT, 1u) + " FT (" + to_string(total_steps) + " steps)");
	print_info("Measurement range: x/h = " + to_string(xh_start, 1u) + " to " + to_string(xh_end, 1u) + " (" + to_string(num_x_stations) + " stations)");
	print_info("Spanwise averaging: ALL " + to_string(z_avg_count) + " z-planes (periodic z, " + to_string(z_avg_count) + "x faster convergence)");
	print_info("----------------------------------------");

	ulong next_report = report_interval;
	bool data_exported = false;
	const float export_interval_FT = 0.5f; // export every 0.5 flow-through after warmup
	float next_export_FT = warmup_FT + export_interval_FT; // first export at 2.5 FT

	// Lambda function to export CSV data
	auto export_csv_data = [&](const string& suffix, bool is_final) {
		const double N = (double)num_total_samples;
		if(N < 1.0) return; // no samples yet

		const string path = get_exe_path();
		print_info("Exporting CSV data" + (suffix.empty() ? "" : " (FT=" + suffix + ")") + " with " + to_string(num_total_samples) + " samples (" + to_string(num_time_samples) + " timesteps x " + to_string(z_avg_count) + " z-planes)...");

		// ===== Figure 7a: Centerline velocity decay (Uc/U0 vs x/h) =====
		{
			const string filename = path + "fig7a_centerline_velocity" + (suffix.empty() ? "" : "_FT" + suffix) + ".csv";
			std::ofstream file(filename);
			file << "x/h,Uc/U0,Uc_rms/U0\n";
			for(int ix = 0; ix < num_x_stations; ix++) {
				if(x_stations[ix] <= 1u) continue; // skip inlet BC region
				if(x_stations[ix] >= Nx - 3u) continue;
				const double Uc_mean = sum_Uc[ix] / N;
				const double Uc_rms = sqrt(fmax(0.0, sum_Uc2[ix] / N - Uc_mean * Uc_mean));
				file << xh_values[ix] << "," << (Uc_mean / lbm_u) << "," << (Uc_rms / lbm_u) << "\n";
			}
			file.close();
			print_info("  Wrote " + filename);
		}

		// ===== Figure 7b: Centerline turbulence intensity (Urms/Uc vs x/h) =====
		{
			const string filename = path + "fig7b_turbulence_intensity" + (suffix.empty() ? "" : "_FT" + suffix) + ".csv";
			std::ofstream file(filename);
			file << "x/h,Urms/Uc\n";
			for(int ix = 0; ix < num_x_stations; ix++) {
				if(x_stations[ix] <= 1u) continue; // skip inlet BC region
				if(x_stations[ix] >= Nx - 3u) continue;
				const double Uc_mean = sum_Uc[ix] / N;
				const double Uc_rms = sqrt(fmax(0.0, sum_Uc2[ix] / N - Uc_mean * Uc_mean));
				const double Urms_over_Uc = (Uc_mean > 1e-10) ? (Uc_rms / Uc_mean) : 0.0;
				file << xh_values[ix] << "," << Urms_over_Uc << "\n";
			}
			file.close();
			print_info("  Wrote " + filename);
		}

		// ===== Half-width (y0.5/h vs x/h) =====
		{
			const string filename = path + "halfwidth" + (suffix.empty() ? "" : "_FT" + suffix) + ".csv";
			std::ofstream file(filename);
			file << "x/h,y0.5/h\n";
			for(int ix = 0; ix < num_x_stations; ix++) {
				if(x_stations[ix] <= 1u) continue; // skip inlet BC region
				if(x_stations[ix] >= Nx - 3u) continue;
				const double Uc_mean = sum_Uc[ix] / N;
				const double U_half = 0.5 * Uc_mean;
				double y_half = 0.0;
				for(uint iy = lateral_extent; iy < lateral_points; iy++) {
					const double U_mean = sum_U_lateral[ix][iy] / N;
					if(U_mean < U_half) {
						const double U_prev = sum_U_lateral[ix][iy-1] / N;
						const double frac = (U_prev - U_half) / (U_prev - U_mean + 1e-10);
						y_half = ((double)(iy - 1 - lateral_extent) + frac);
						break;
					}
				}
				file << xh_values[ix] << "," << (y_half / (double)h_cells) << "\n";
			}
			file.close();
			print_info("  Wrote " + filename);
		}

		// ===== Figure 4: Lateral profiles (U/Uc vs y/y0.5) =====
		// Ahmed et al. uses y/y0.5 (normalized by local half-width) for self-similar profiles
		{
			// x/h stations: 0, 2, 4, 6, 8, 10
			const float fig4_xh[] = {0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f};
			const int num_fig4 = 6;
			int fig4_stations[6];
			for(int i = 0; i < num_fig4; i++) {
				fig4_stations[i] = (int)round((fig4_xh[i] - xh_start) / xh_increment);
			}
			const string filename = path + "fig4_lateral_profiles" + (suffix.empty() ? "" : "_FT" + suffix) + ".csv";
			std::ofstream file(filename);
			file << "y/y0.5"; // Changed from y/h to y/y0.5 to match Ahmed et al. Fig 4
			for(int i = 0; i < num_fig4; i++) {
				file << ",U/Uc_xh" << (int)fig4_xh[i];
			}
			file << "\n";
			// Compute Uc and y0.5 at each station
			vector<double> Uc_at_station(num_fig4);
			vector<double> y05_at_station(num_fig4);
			for(int i = 0; i < num_fig4; i++) {
				int ix = fig4_stations[i];
				if(ix >= 0 && ix < num_x_stations && x_stations[ix] > 1u && x_stations[ix] < Nx - 3u) {
					Uc_at_station[i] = sum_Uc[ix] / N;
					// Find local half-width for this station
					const double U_half = 0.5 * Uc_at_station[i];
					y05_at_station[i] = 0.5 * (double)h_cells; // default: 0.5h
					for(uint iy2 = lateral_extent; iy2 < lateral_points; iy2++) {
						const double U_mean2 = sum_U_lateral[ix][iy2] / N;
						if(U_mean2 < U_half) {
							const double U_prev2 = sum_U_lateral[ix][iy2-1] / N;
							const double frac2 = (U_prev2 - U_half) / (U_prev2 - U_mean2 + 1e-10);
							y05_at_station[i] = (double)(iy2 - 1 - lateral_extent) + frac2;
							break;
						}
					}
					if(y05_at_station[i] < 0.1) y05_at_station[i] = 0.5 * (double)h_cells; // safety floor
				} else {
					Uc_at_station[i] = lbm_u;
					y05_at_station[i] = 0.5 * (double)h_cells;
				}
			}
			for(uint iy = 0; iy < lateral_points; iy++) {
				const double y_rel = (double)iy - (double)lateral_extent;
				// Station-specific y/y0.5 normalization for proper self-similar collapse
				// Each station uses its OWN y0.5, matching Ahmed et al. Fig 4 convention
				const double y_over_y05 = y_rel / y05_at_station[0]; // first column: y/y0.5 from x/h=0 as reference
				file << y_over_y05;
				for(int i = 0; i < num_fig4; i++) {
					int ix = fig4_stations[i];
					if(ix >= 0 && ix < num_x_stations && x_stations[ix] > 1u && x_stations[ix] < Nx - 3u) {
						const double U_mean = sum_U_lateral[ix][iy] / N;
						// U/Uc with station-specific normalization
						file << "," << (U_mean / Uc_at_station[i]);
					} else {
						file << ",";
					}
				}
				file << "\n";
			}
			// Also write station-specific y0.5/h values for post-processing
			// Analysis scripts can compute station-specific y/y0.5 = (y/h) / (y0.5/h)
			file.close();
			{
				const string y05_filename = path + "fig4_y05_values" + (suffix.empty() ? "" : "_FT" + suffix) + ".csv";
				std::ofstream y05_file(y05_filename);
				y05_file << "x/h,y0.5/h,y0.5_cells\n";
				for(int i = 0; i < num_fig4; i++) {
					y05_file << fig4_xh[i] << "," << (y05_at_station[i] / (double)h_cells) << "," << y05_at_station[i] << "\n";
				}
				y05_file.close();
				print_info("  Wrote " + y05_filename);
			}
			print_info("  Wrote " + filename);
		}

		// ===== Figure 10: Volume flow rate (Q*/Q0* vs x/h) =====
		// Per Ahmed et al.: Q* = Q / y0.5 (specific volume flow rate normalized by half-width)
		// This removes the jet spreading effect from the entrainment measurement
		{
			// First, calculate half-width y0.5 for each station (same logic as halfwidth export)
			vector<double> y_half_values(num_x_stations, 0.0);
			for(int ix = 0; ix < num_x_stations; ix++) {
				if(x_stations[ix] <= 1u) continue;
				if(x_stations[ix] >= Nx - 3u) continue;
				const double Uc_mean = sum_Uc[ix] / N;
				const double U_half = 0.5 * Uc_mean;
				for(uint iy = lateral_extent; iy < lateral_points; iy++) {
					const double U_mean = sum_U_lateral[ix][iy] / N;
					if(U_mean < U_half) {
						const double U_prev = sum_U_lateral[ix][iy-1] / N;
						const double frac = (U_prev - U_half) / (U_prev - U_mean + 1e-10);
						y_half_values[ix] = ((double)(iy - 1 - lateral_extent) + frac);
						break;
					}
				}
				// Ensure minimum half-width of 0.5*h to avoid division issues
				if(y_half_values[ix] < 0.5 * (double)h_cells) {
					y_half_values[ix] = 0.5 * (double)h_cells;
				}
			}

			// Find Q0* at first valid station (x/h >= 0)
			// Q* = Q / y0.5 (specific volume flow rate)
			const double Q_threshold = 0.1 * lbm_u * (double)h_cells;
			double Q0_star = 0.0;
			for(int ix = 0; ix < num_x_stations; ix++) {
				if(x_stations[ix] <= 1u) continue;
				if(x_stations[ix] >= Nx - 3u) continue;
				if(xh_values[ix] < 0.0f) continue; // Q0* should be at x/h >= 0
				double Q_mean = sum_Q[ix] / N;
				double y_half = y_half_values[ix];
				if(Q_mean > Q_threshold && y_half > 0.0) {
					Q0_star = Q_mean / y_half; // Q* = Q / y0.5
					break;
				}
			}

			// If no valid Q0* found, skip output with warning
			if(Q0_star <= 0.0) {
				print_info("  WARNING: No valid Q0* reference found, skipping volume flow output");
			} else {
				const string filename = path + "fig10_volume_flow" + (suffix.empty() ? "" : "_FT" + suffix) + ".csv";
				std::ofstream file(filename);
				file << "x/h,Qstar/Q0star\n"; // Changed header to reflect normalized quantity
				for(int ix = 0; ix < num_x_stations; ix++) {
					if(x_stations[ix] <= 1u) continue;
					if(x_stations[ix] >= Nx - 3u) continue;
					const double Q_mean = sum_Q[ix] / N;
					const double y_half = y_half_values[ix];
					if(Q_mean > 0.0 && y_half > 0.0) {
						const double Q_star = Q_mean / y_half; // Q* = Q / y0.5
						file << xh_values[ix] << "," << (Q_star / Q0_star) << "\n";
					}
				}
				file.close();
				print_info("  Wrote " + filename);
			}
		}

		if(is_final) {
			print_info("Final data export complete.");
		}
	};

	while(lbm.get_t() <= total_steps) {
		lbm.run(sample_interval);

		const float ft = (float)lbm.get_t() / flow_through_time;
		const bool in_averaging_phase = (lbm.get_t() > warmup_steps);

		// ============ Collect statistics during averaging phase ============
		// Uses spanwise averaging: sample at ALL Nz z-planes (periodic z)
		if(in_averaging_phase) {
			lbm.u.read_from_device();
			num_time_samples++;
			num_total_samples += (ulong)Nz; // all z-planes

			for(int ix = 0; ix < num_x_stations; ix++) {
				const uint x = x_stations[ix];
				if(x <= 1u) continue; // skip inlet BC region (x=0 is TYPE_X inlet)
				if(x >= Nx - 3u) continue; // skip outlet region

				// Loop over ALL z-planes (periodic z → all planes are valid)
				for(uint z = 0u; z < Nz; z++) {

					// Centerline velocity (Fig 7a) - at y=cy for each z-plane
					const ulong n_center = (ulong)x + ((ulong)cy + (ulong)z * (ulong)Ny) * (ulong)Nx;
					const float Uc = lbm.u.x[n_center];
					sum_Uc[ix] += (double)Uc;
					sum_Uc2[ix] += (double)Uc * (double)Uc; // double precision squaring for variance accuracy

					// Lateral profile (Fig 4) - at each y for this z-plane
					for(uint iy = 0; iy < lateral_points; iy++) {
						const uint y = cy - lateral_extent + iy;
						if(y >= Ny) continue;
						const ulong n = (ulong)x + ((ulong)y + (ulong)z * (ulong)Ny) * (ulong)Nx;
						const float U = lbm.u.x[n];
						sum_U_lateral[ix][iy] += (double)U;
						sum_U2_lateral[ix][iy] += (double)U * (double)U; // double precision squaring
					}

					// Volume flow rate (Fig 10) - 2D integration over jet cross-section (y,z)
					// Per Ahmed et al.: Q = integral(U dy dz) with 10% velocity cutoff
					// Each z-plane contributes one z-slice to the full 2D integral
					const float Uc_local = lbm.u.x[n_center];
					const float cutoff = 0.1f * fabs(Uc_local); // 10% of centerline velocity

					// Skip stations where centerline velocity is negligible
					if(fabs(Uc_local) >= 1e-6f) {
						double Q_slice = 0.0; // Q contribution from this z-plane (1D y-integral)
						// Integrate from centerline outward in +y direction until U < cutoff
						for(uint y = cy; y < Ny; y++) {
							const ulong n = (ulong)x + ((ulong)y + (ulong)z * (ulong)Ny) * (ulong)Nx;
							const float U = lbm.u.x[n];
							if(fabs(U) < cutoff) break;
							Q_slice += (double)U;
						}
						// Integrate from centerline outward in -y direction until U < cutoff
						for(int y = (int)cy - 1; y >= 0; y--) {
							const ulong n = (ulong)x + ((ulong)y + (ulong)z * (ulong)Ny) * (ulong)Nx;
							const float U = lbm.u.x[n];
							if(fabs(U) < cutoff) break;
							Q_slice += (double)U;
						}
						sum_Q[ix] += Q_slice; // accumulate z-slices for 2D integral
					}
				}
			}
		}

		// ============ Periodic reporting ============
		if(lbm.get_t() >= next_report) {
			next_report += report_interval;

			lbm.u.read_from_device();

			// Quick centerline check (only show x/h >= 0 for readability)
			string status = in_averaging_phase ? "AVERAGING" : "WARMUP";
			string centerline_str = "";
			for(int ix = 0; ix < num_x_stations; ix += 5) { // every 5th station
				if(x_stations[ix] <= 1u) continue; // skip inlet BC region
				if(x_stations[ix] >= Nx - 3u) continue;
				const ulong n = (ulong)x_stations[ix] + ((ulong)cy + (ulong)cz * (ulong)Ny) * (ulong)Nx;
				const float Uc = lbm.u.x[n];
				centerline_str += " x/h=" + to_string((int)xh_values[ix]) + ":" + to_string(Uc/lbm_u, 2u);
			}

			print_info("FT=" + to_string(ft, 2u) + " [" + status + "] samples=" + to_string(num_time_samples) + "x" + to_string(z_avg_count) + "=" + to_string(num_total_samples) + " | Uc/U0:" + centerline_str);

#ifdef VREMAN_SGS
			// Monitor Vreman nu_sgs at key centerline locations
			string vreman_str = "";
			const int monitor_xh[] = {0, 2, 5, 10, 15, 20}; // x/h locations to monitor
			for(int i = 0; i < 6; i++) {
				const int xh = monitor_xh[i];
				const uint x_pos = (uint)(xh * h_cells);
				if(x_pos >= 3 && x_pos < Nx - 3) {
					const float nu_sgs = compute_vreman_nu_sgs(lbm, x_pos, cy, cz);
					vreman_str += " x/h=" + to_string(xh) + ":" + to_string(nu_sgs, 6u);
				}
			}
			print_info("Vreman nu_sgs:" + vreman_str);
#endif // VREMAN_SGS
		}

		// ============ Periodic intermediate export during averaging phase ============
		if(in_averaging_phase && ft >= next_export_FT && !data_exported) {
			print_info("----------------------------------------");
			export_csv_data(to_string(next_export_FT, 1u), false);
			print_info("----------------------------------------");
			next_export_FT += export_interval_FT;
		}

		// ============ Final export after averaging phase ============
		if(lbm.get_t() >= total_steps && !data_exported) {
			print_info("========================================");
			print_info("FINAL EXPORT at FT=" + to_string(ft, 2u));
			export_csv_data("", true); // empty suffix for final files
			print_info("========================================");
			data_exported = true;
		}
	}

	// Stop simulation after export is complete
	print_info("========================================");
	print_info("Simulation complete at FT=35. Exiting...");
	print_info("========================================");
	// lbm.run(); // Commented out: auto-stop after 35 FT instead of running indefinitely
}
#endif // !BENCHMARK
/**/



/*void main_setup() { // 3D Taylor-Green vortices; required extensions in defines.hpp: INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(128u, 128u, 128u, 1u, 1u, 1u, 0.01f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		const float A = 0.25f;
		const uint periodicity = 1u;
		const float a=(float)Nx/(float)periodicity, b=(float)Ny/(float)periodicity, c=(float)Nz/(float)periodicity;
		const float fx = (float)x+0.5f-0.5f*(float)Nx;
		const float fy = (float)y+0.5f-0.5f*(float)Ny;
		const float fz = (float)z+0.5f-0.5f*(float)Nz;
		lbm.u.x[n] =  A*cosf(2.0f*pif*fx/a)*sinf(2.0f*pif*fy/b)*sinf(2.0f*pif*fz/c);
		lbm.u.y[n] = -A*sinf(2.0f*pif*fx/a)*cosf(2.0f*pif*fy/b)*sinf(2.0f*pif*fz/c);
		lbm.u.z[n] =  A*sinf(2.0f*pif*fx/a)*sinf(2.0f*pif*fy/b)*cosf(2.0f*pif*fz/c);
		lbm.rho[n] = 1.0f-sq(A)*3.0f/4.0f*(cosf(4.0f*pif*fx/a)+cosf(4.0f*pif*fy/b));
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_STREAMLINES;
	lbm.run();
	//lbm.run(1000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/2u, Ny/2u, Nz/2u)]); wait(); // test for binary identity
} /**/



/*void main_setup() { // 2D Taylor-Green vortices (use D2Q9); required extensions in defines.hpp: INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(1024u, 1024u, 1u, 0.02f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		const float A = 0.2f;
		const uint periodicity = 5u;
		const float a=(float)Nx/(float)periodicity, b=(float)Ny/(float)periodicity;
		const float fx = (float)x+0.5f-0.5f*(float)Nx;
		const float fy = (float)y+0.5f-0.5f*(float)Ny;
		lbm.u.x[n] =  A*cosf(2.0f*pif*fx/a)*sinf(2.0f*pif*fy/b);
		lbm.u.y[n] = -A*sinf(2.0f*pif*fx/a)*cosf(2.0f*pif*fy/b);
		lbm.rho[n] = 1.0f-sq(A)*3.0f/4.0f*(cosf(4.0f*pif*fx/a)+cosf(4.0f*pif*fy/b));
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FIELD;
	lbm.graphics.slice_mode = 3;
	lbm.run();
} /**/



/*void main_setup() { // Poiseuille flow validation; required extensions in defines.hpp: VOLUME_FORCE
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint R = 63u; // channel radius (default: 63)
	const float umax = 0.1f; // maximum velocity in channel center (must be < 0.57735027f)
	const float tau = 1.0f; // relaxation time (must be > 0.5f), tau = nu*3+0.5
	const float nu = units.nu_from_tau(tau); // nu = (tau-0.5)/3
	const uint H = 2u*(R+1u);
#ifndef D2Q9
	LBM lbm(H, lcm(sq(H), WORKGROUP_SIZE)/sq(H), H, nu, 0.0f, units.f_from_u_Poiseuille_3D(umax, 1.0f, nu, R), 0.0f); // 3D
#else // D2Q9
	LBM lbm(lcm(H, WORKGROUP_SIZE)/H, H, 1u, nu, units.f_from_u_Poiseuille_2D(umax, 1.0f, nu, R), 0.0f, 0.0f); // 2D
#endif // D2Q9
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
#ifndef D2Q9
		if(!cylinder(x, y, z, lbm.center(), float3(0u, Ny, 0u), 0.5f*(float)min(Nx, Nz)-1.0f)) lbm.flags[n] = TYPE_S; // 3D
#else // D2Q9
		if(y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S; // 2D
#endif // D2Q9
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	double error_min = max_double;
	while(true) { // main simulation loop
		lbm.run(1000u);
		lbm.u.read_from_device();
		double error_dif=0.0, error_sum=0.0;
#ifndef D2Q9
		for(uint x=0u; x<Nx; x++) {
			for(uint y=Ny/2u; y<Ny/2u+1u; y++) {
				for(uint z=0; z<Nz; z++) {
					const uint n = x+(y+z*Ny)*Nx;
					const double r = (double)sqrt(sq(x+0.5f-0.5f*(float)Nx)+sq(z+0.5f-0.5f*(float)Nz)); // radius from channel center
					if(r<R) {
						const double unum = (double)sqrt(sq(lbm.u.x[n])+sq(lbm.u.y[n])+sq(lbm.u.z[n])); // numerical velocity
						const double uref = umax*(sq(R)-sq(r))/sq(R); // theoretical velocity profile u = G*(R^2-r^2)
						error_dif += sq(unum-uref); // L2 error (Krüger p. 138)
						error_sum += sq(uref);
					}
				}
			}
		}
#else // D2Q9
		for(uint x=Nx/2u; x<Nx/2u+1u; x++) {
			for(uint y=1u; y<Ny-1u; y++) {
				const uint n = x+(y+0u*Ny)*Nx;
				const double r = (double)(y+0.5f-0.5f*(float)Ny); // radius from channel center
				const double unum = (double)sqrt(sq(lbm.u.x[n])+sq(lbm.u.y[n])); // numerical velocity
				const double uref = umax*(sq(R)-sq(r))/sq(R); // theoretical velocity profile u = G*(R^2-r^2)
				error_dif += sq(unum-uref); // L2 error (Krüger p. 138)
				error_sum += sq(uref);
			}
		}
#endif // D2Q9
		if(sqrt(error_dif/error_sum)>=error_min) { // stop when error has converged
			print_info("Poiseuille flow error converged after "+to_string(lbm.get_t())+" steps to "+to_string(100.0*error_min, 3u)+"%"); // typical expected L2 errors: 2-5% (Krüger p. 256)
			wait();
			exit(0);
		}
		error_min = fmin(error_min, sqrt(error_dif/error_sum));
		print_info("Poiseuille flow error after t="+to_string(lbm.get_t())+" is "+to_string(100.0*error_min, 3u)+"%"); // typical expected L2 errors: 2-5% (Krüger p. 256)
	}
} /**/



/*void main_setup() { // Stokes drag validation; required extensions in defines.hpp: FORCE_FIELD, EQUILIBRIUM_BOUNDARIES
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const ulong dt = 100ull; // check error every dt time steps
	const float R = 32.0f; // sphere radius
	const float Re = 0.01f; // Reynolds number
	const float nu = 1.0f; // kinematic shear viscosity
	const float rho = 1.0f; // density
	const uint L = to_uint(8.0f*R); // simulation box size
	const float u = units.u_from_Re(Re, 2.0f*R, nu); // velocity
	LBM lbm(L, L, L, nu); // flow driven by equilibrium boundaries
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E;
		if(sphere(x, y, z, lbm.center(), R)) {
			lbm.flags[n] = TYPE_S|TYPE_X; // flag boundary cells for force summation additionally with TYPE_X
		} else {
			lbm.rho[n] = units.rho_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R, rho, nu);
			const float3 un = units.u_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R);
			lbm.u.x[n] = un.x;
			lbm.u.y[n] = un.y;
			lbm.u.z[n] = un.z;
		}
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	double E1=1000.0, E2=1000.0;
	while(true) { // main simulation loop
		lbm.run(dt);
		const float3 force = lbm.object_force(TYPE_S|TYPE_X);
		const double F_theo = units.F_Stokes(rho, u, nu, R);
		const double F_sim = (double)length(force);
		const double E0 = fabs(F_sim-F_theo)/F_theo;
		print_info(to_string(lbm.get_t())+", expected: "+to_string(F_theo, 6u)+", measured: "+to_string(F_sim, 6u)+", error = "+to_string((float)(100.0*E0), 1u)+"%");
		if(converged(E2, E1, E0, 1E-4)) { // stop when error has sufficiently converged
			print_info("Error converged after "+to_string(lbm.get_t())+" steps to "+to_string(100.0*E0, 1u)+"%");
			wait();
			break;
		}
		E2 = E1;
		E1 = E0;
	}
} /**/



/*void main_setup() { // cylinder in rectangular duct; required extensions in defines.hpp: VOLUME_FORCE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const float Re = 25000.0f;
	const float D = 64.0f;
	const float u = rsqrt(3.0f);
	const float w=D, l=12.0f*D, h=3.0f*D;
	const float nu = units.nu_from_Re(Re, D, u);
	const float f = units.f_from_u_rectangular_duct(w, D, 1.0f, nu, u);
	LBM lbm(to_uint(w), to_uint(l), to_uint(h), nu, 0.0f, f, 0.0f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		lbm.u.y[n] = 0.1f*u;
		if(cylinder(x, y, z, float3(lbm.center().x, 2.0f*D, lbm.center().z), float3(Nx, 0u, 0u), 0.5f*D)) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // x and z non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_Q_CRITERION;
	lbm.run();
} /**/



/*void main_setup() { // Taylor-Couette flow; required extensions in defines.hpp: MOVING_BOUNDARIES, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(96u, 96u, 192u, 1u, 1u, 1u, 0.04f);
	// ###################################################################################### define geometry ######################################################################################
	const uint threads = (uint)thread::hardware_concurrency();
	vector<uint> seed(threads);
	for(uint t=0u; t<threads; t++) seed[t] = 42u+t;
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), threads, [&](ulong n, uint t) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(!cylinder(x, y, z, lbm.center(), float3(0u, 0u, Nz), (float)(Nx/2u-1u))) lbm.flags[n] = TYPE_S;
		if( cylinder(x, y, z, lbm.center(), float3(0u, 0u, Nz), (float)(Nx/4u   ))) {
			const float3 relative_position = lbm.relative_position(n);
			lbm.u.x[n] =  relative_position.y;
			lbm.u.y[n] = -relative_position.x;
			lbm.u.z[n] = (1.0f-random(seed[t], 2.0f))*0.001f;
			lbm.flags[n] = TYPE_S;
		}
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_STREAMLINES;
	lbm.run();
	//lbm.run(4000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/4u, Ny/4u, Nz/2u)]); wait(); // test for binary identity
} /**/



/*void main_setup() { // lid-driven cavity; required extensions in defines.hpp: MOVING_BOUNDARIES, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint L = 128u;
	const float Re = 1000.0f;
	const float u = 0.1f;
	LBM lbm(L, L, L, units.nu_from_Re(Re, (float)(L-2u), u));
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z==Nz-1) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_STREAMLINES;
	lbm.run();
} /**/



/*void main_setup() { // 2D Karman vortex street; required extensions in defines.hpp: D2Q9, FP16S, EQUILIBRIUM_BOUNDARIES, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint R = 16u;
	const float Re = 250.0f;
	const float u = 0.10f;
	LBM lbm(16u*R, 32u*R, 1u, units.nu_from_Re(Re, 2.0f*(float)R, u));
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(cylinder(x, y, z, float3(Nx/2u, Ny/4u, Nz/2u), float3(0u, 0u, Nz), (float)R)) lbm.flags[n] = TYPE_S;
		else lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_FIELD;
	lbm.graphics.slice_mode = 3;
	lbm.run();
} /**/



/*void main_setup() { // particle test; required extensions in defines.hpp: VOLUME_FORCE, FORCE_FIELD, MOVING_BOUNDARIES, PARTICLES, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint L = 128u;
	const float Re = 1000.0f;
	const float u = 0.1f;
	LBM lbm(L, L, L, units.nu_from_Re(Re, (float)(L-2u), u), 0.0f, 0.0f, -0.00001f, cb(L/4u), 2.0f);
	// ###################################################################################### define geometry ######################################################################################
	uint seed = 42u;
	for(ulong n=0ull; n<lbm.particles->length(); n++) {
		lbm.particles->x[n] = random_symmetric(seed, 0.5f*lbm.size().x/4.0f);
		lbm.particles->y[n] = random_symmetric(seed, 0.5f*lbm.size().y/4.0f);
		lbm.particles->z[n] = random_symmetric(seed, 0.5f*lbm.size().z/4.0f);
	}
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z==Nz-1) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_STREAMLINES|VIS_PARTICLES;
	lbm.run();
} /**/



/*void main_setup() { // delta wing; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint L = 128u;
	const float Re = 100000.0f;
	const float u = 0.075f;
	LBM lbm(L, 4u*L, L, units.nu_from_Re(Re, (float)L, u));
	// ###################################################################################### define geometry ######################################################################################
	const float3 offset = float3(lbm.center().x, 0.0f, lbm.center().z);
	const float3 p0 = offset+float3(  0*(int)L/64,  5*(int)L/64,  20*(int)L/64);
	const float3 p1 = offset+float3(-20*(int)L/64, 90*(int)L/64, -10*(int)L/64);
	const float3 p2 = offset+float3(+20*(int)L/64, 90*(int)L/64, -10*(int)L/64);
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(triangle(x, y, z, p0, p1, p2)) lbm.flags[n] = TYPE_S;
		else lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run();
} /**/



/*void main_setup() { // NASA Common Research Model; required extensions in defines.hpp: FP16C, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 1.5f, 1.0f/3.0f), 2000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float Re = 10000000.0f;
	const float u = 0.075f;
	LBM lbm(lbm_N, units.nu_from_Re(Re, (float)lbm_N.x, u));
	// ###################################################################################### define geometry ######################################################################################
	// model: https://commonresearchmodel.larc.nasa.gov/high-lift-crm/high-lift-crm-geometry/assembled-geometry/, .stp file converted to .stl with https://imagetostl.com/convert/file/stp/to/stl
	Mesh* half = read_stl(get_exe_path()+"../stl/crm-hl_reference_ldg.stl", lbm.size(), float3(0.0f), float3x3(float3(0, 0, 1), radians(90.0f)), 1.0f*lbm_N.x);
	half->translate(float3(-0.5f*(half->pmax.x-half->pmin.x), 0.0f, 0.0f));
	Mesh* mesh = new Mesh(2u*half->triangle_number, float3(0.0f));
	for(uint i=0u; i<half->triangle_number; i++) {
		mesh->p0[i] = half->p0[i];
		mesh->p1[i] = half->p1[i];
		mesh->p2[i] = half->p2[i];
	}
	half->rotate(float3x3(float3(1, 0, 0), radians(180.0f))); // mirror-copy half
	for(uint i=0u; i<half->triangle_number; i++) {
		mesh->p0[half->triangle_number+i] = -half->p0[i];
		mesh->p1[half->triangle_number+i] = -half->p1[i];
		mesh->p2[half->triangle_number+i] = -half->p2[i];
	}
	delete half;
	mesh->find_bounds();
	mesh->rotate(float3x3(float3(1, 0, 0), radians(-10.0f)));
	mesh->translate(float3(0.0f, 0.0f, -0.5f*(mesh->pmin.z+mesh->pmax.z)));
	mesh->translate(float3(0.0f, -0.5f*lbm.size().y+mesh->pmax.y+0.5f*(lbm.size().x-(mesh->pmax.x-mesh->pmin.x)), 0.0f));
	mesh->translate(lbm.center());
	lbm.voxelize_mesh_on_device(mesh);
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run();
} /**/



/*void main_setup() { // Concorde; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 3.0f, 0.5f), 2084u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float si_u = 300.0f/3.6f;
	const float si_length=62.0f, si_width=26.0f;
	const float si_T = 1.0f;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	const float lbm_length = 0.56f*(float)lbm_N.y;
	const float lbm_u = 0.075f;
	units.set_m_kg_s(lbm_length, lbm_u, 1.0f, si_length, si_u, si_rho);
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	print_info("Re = "+to_string(to_uint(units.si_Re(si_width, si_u, si_nu))));
	LBM lbm(lbm_N, 1u, 1u, 1u, lbm_nu);
	// ###################################################################################### define geometry ######################################################################################
	const float3 center = float3(lbm.center().x, 0.52f*lbm_length, lbm.center().z+0.03f*lbm_length);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-10.0f))*float3x3(float3(0, 0, 1), radians(90.0f))*float3x3(float3(1, 0, 0), radians(90.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/concord_cut_large.stl", center, rotation, lbm_length); // https://www.thingiverse.com/thing:1176931/files
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run(0u, lbm_T); // initialize simulation
	lbm.write_status();
	while(lbm.get_t()<=lbm_T) { // main simulation loop
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		if(lbm.graphics.next_frame(lbm_T, 10.0f)) {
			lbm.graphics.set_camera_free(float3(0.491343f*(float)Nx, -0.882147f*(float)Ny, 0.564339f*(float)Nz), -78.0f, 6.0f, 22.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/front/");
			lbm.graphics.set_camera_free(float3(1.133361f*(float)Nx, 1.407077f*(float)Ny, 1.684411f*(float)Nz), 72.0f, 12.0f, 20.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/back/");
			lbm.graphics.set_camera_centered(0.0f, 0.0f, 25.0f, 1.648722f);
			lbm.graphics.write_frame(get_exe_path()+"export/side/");
			lbm.graphics.set_camera_centered(0.0f, 90.0f, 25.0f, 1.648722f);
			lbm.graphics.write_frame(get_exe_path()+"export/top/");
			lbm.graphics.set_camera_free(float3(0.269361f*(float)Nx, -0.179720f*(float)Ny, 0.304988f*(float)Nz), -56.0f, 31.6f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/wing/");
			lbm.graphics.set_camera_free(float3(0.204399f*(float)Nx, 0.340055f*(float)Ny, 1.620902f*(float)Nz), 80.0f, 35.6f, 34.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/follow/");
		}
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
		lbm.run(1u, lbm_T); // run dt time steps
	}
	lbm.write_status();
} /**/



/*void main_setup() { // Boeing 747; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 0.5f), 880u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_Re = 1000000.0f;
	const float lbm_u = 0.075f;
	const ulong lbm_T = 10000ull;
	LBM lbm(lbm_N, units.nu_from_Re(lbm_Re, (float)lbm_N.x, lbm_u));
	// ###################################################################################### define geometry ######################################################################################
	const float size = 1.0f*lbm.size().x;
	const float3 center = float3(lbm.center().x, 0.55f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-15.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/techtris_airplane.stl", center, rotation, size); // https://www.thingiverse.com/thing:2772812/files
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
	lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 10.0f)) lbm.graphics.write_frame(); // render enough frames 10 seconds of 60fps video
		lbm.run(1u, lbm_T);
	}
#else // GRAPHICS && !INTERACTIVE_GRAPHICS
	lbm.run();
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
} /**/



/*void main_setup() { // Star Wars X-wing; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 0.5f), 880u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_Re = 100000.0f;
	const float lbm_u = 0.075f;
	const ulong lbm_T = 50000ull;
	LBM lbm(lbm_N, units.nu_from_Re(lbm_Re, (float)lbm_N.x, lbm_u));
	// ###################################################################################### define geometry ######################################################################################
	const float size = 1.0f*lbm.size().x;
	const float3 center = float3(lbm.center().x, 0.55f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(0, 0, 1), radians(180.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/X-Wing.stl", center, rotation, size); // https://www.thingiverse.com/thing:353276/files
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 30.0f)) {
			lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/t/");
			lbm.graphics.set_camera_free(float3(0.5f*(float)Nx, -0.35f*(float)Ny, -0.7f*(float)Nz), -33.0f, -40.0f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/b/");
			lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.51f*(float)Ny, 0.75f*(float)Nz), 90.0f, 28.0f, 80.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/f/");
			lbm.graphics.set_camera_free(float3(0.7f*(float)Nx, -0.15f*(float)Ny, 0.06f*(float)Nz), 0.0f, 0.0f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/s/");
		}
		lbm.run(1u, lbm_T);
	}
#else // GRAPHICS && !INTERACTIVE_GRAPHICS
	lbm.run();
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
} /**/



/*void main_setup() { // Star Wars TIE fighter; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 1.0f), 1760u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_Re = 100000.0f;
	const float lbm_u = 0.075f;
	const ulong lbm_T = 50000ull;
	const ulong lbm_dt = 28ull;
	LBM lbm(lbm_N, units.nu_from_Re(lbm_Re, (float)lbm_N.x, lbm_u));
	// ###################################################################################### define geometry ######################################################################################
	const float size = 0.65f*lbm.size().x;
	const float3 center = float3(lbm.center().x, 0.6f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(90.0f));
	Mesh* mesh = read_stl(get_exe_path()+"../stl/DWG_Tie_Fighter_Assembled_02.stl", lbm.size(), center, rotation, size); // https://www.thingiverse.com/thing:2919109/files
	lbm.voxelize_mesh_on_device(mesh);
	lbm.flags.read_from_device();
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<lbm_T) { // main simulation loop
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		if(lbm.graphics.next_frame(lbm_T, 30.0f)) {
			lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 0.63f*(float)Nz), -33.0f, 33.0f, 80.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/t/");
			lbm.graphics.set_camera_free(float3(0.3f*(float)Nx, -1.5f*(float)Ny, -0.45f*(float)Nz), -83.0f, -10.0f, 40.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/b/");
			lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.57f*(float)Ny, 0.7f*(float)Nz), 90.0f, 29.5f, 80.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/f/");
			lbm.graphics.set_camera_free(float3(2.5f*(float)Nx, 0.0f*(float)Ny, 0.0f*(float)Nz), 0.0f, 0.0f, 50.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/s/");
		}
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
		lbm.run(lbm_dt, lbm_T);
		const float3x3 rotation = float3x3(float3(0.2f, 1.0f, 0.1f), radians(0.4032f)); // create rotation matrix to rotate mesh
		lbm.unvoxelize_mesh_on_device(mesh);
		mesh->rotate(rotation); // rotate mesh
		lbm.voxelize_mesh_on_device(mesh);
	}
} /**/



/*void main_setup() { // radial fan; required extensions in defines.hpp: FP16S, MOVING_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(3.0f, 3.0f, 1.0f), 181u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_Re = 100000.0f;
	const float lbm_u = 0.1f;
	const ulong lbm_T = 48000ull;
	const ulong lbm_dt = 10ull;
	LBM lbm(lbm_N, 1u, 1u, 1u, units.nu_from_Re(lbm_Re, (float)lbm_N.x, lbm_u));
	// ###################################################################################### define geometry ######################################################################################
	const float radius = 0.25f*(float)lbm_N.x;
	const float3 center = float3(lbm.center().x, lbm.center().y, 0.36f*radius);
	const float lbm_omega=lbm_u/radius, lbm_domega=lbm_omega*lbm_dt;
	Mesh* mesh = read_stl(get_exe_path()+"../stl/FAN_Solid_Bottom.stl", lbm.size(), center, 2.0f*radius); // https://www.thingiverse.com/thing:6113/files
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u) lbm.flags[n] = TYPE_S; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<lbm_T) { // main simulation loop
		lbm.voxelize_mesh_on_device(mesh, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, lbm_omega));
		lbm.run(lbm_dt, lbm_T);
		mesh->rotate(float3x3(float3(0.0f, 0.0f, 1.0f), lbm_domega)); // rotate mesh
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		if(lbm.graphics.next_frame(lbm_T, 30.0f)) {
			lbm.graphics.set_camera_free(float3(0.353512f*(float)Nx, -0.150326f*(float)Ny, 1.643939f*(float)Nz), -25.0f, 61.0f, 100.0f);
			lbm.graphics.write_frame();
		}
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
	}
} /**/



/*void main_setup() { // electric ducted fan (EDF); required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, MOVING_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 1.5f, 1.0f), 8000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_Re = 1000000.0f;
	const float lbm_u = 0.1f;
	const ulong lbm_T = 180000ull;
	const ulong lbm_dt = 4ull;
	LBM lbm(lbm_N, units.nu_from_Re(lbm_Re, (float)lbm_N.x, lbm_u));
	// ###################################################################################### define geometry ######################################################################################
	const float3 center = lbm.center();
	const float3x3 rotation = float3x3(float3(0, 0, 1), radians(180.0f));
	Mesh* stator = read_stl(get_exe_path()+"../stl/edf_v39.stl", 1.0f, rotation); // https://www.thingiverse.com/thing:3014759/files
	Mesh* rotor = read_stl(get_exe_path()+"../stl/edf_v391.stl", 1.0f, rotation); // https://www.thingiverse.com/thing:3014759/files
	const float scale = 0.98f*stator->get_scale_for_box_fit(lbm.size()); // scale stator and rotor to simulation box size
	stator->scale(scale);
	rotor->scale(scale);
	stator->translate(lbm.center()-stator->get_bounding_box_center()-float3(0.0f, 0.2f*stator->get_max_size(), 0.0f)); // move stator and rotor to simulation box center
	rotor->translate(lbm.center()-rotor->get_bounding_box_center()-float3(0.0f, 0.41f*stator->get_max_size(), 0.0f));
	stator->set_center(stator->get_center_of_mass()); // set rotation center of mesh to its center of mass
	rotor->set_center(rotor->get_center_of_mass());
	const float lbm_radius=0.5f*rotor->get_max_size(), omega=lbm_u/lbm_radius, domega=omega*(float)lbm_dt;
	lbm.voxelize_mesh_on_device(stator, TYPE_S, center);
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]==0u) lbm.u.y[n] = 0.3f*lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<lbm_T) { // main simulation loop
		lbm.voxelize_mesh_on_device(rotor, TYPE_S, center, float3(0.0f), float3(0.0f, omega, 0.0f));
		lbm.run(lbm_dt, lbm_T);
		rotor->rotate(float3x3(float3(0.0f, 1.0f, 0.0f), domega)); // rotate mesh
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		if(lbm.graphics.next_frame(lbm_T, 30.0f)) {
			lbm.graphics.set_camera_centered(-70.0f+100.0f*(float)lbm.get_t()/(float)lbm_T, 2.0f, 60.0f, 1.284025f);
			lbm.graphics.write_frame();
		}
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
	}
} /**/



/*void main_setup() { // aerodynamics of a cow; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 1.0f), 1000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float si_u = 1.0f;
	const float si_length = 2.4f;
	const float si_T = 10.0f;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	const float lbm_length = 0.65f*(float)lbm_N.y;
	const float lbm_u = 0.075f;
	units.set_m_kg_s(lbm_length, lbm_u, 1.0f, si_length, si_u, si_rho);
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	print_info("Re = "+to_string(to_uint(units.si_Re(si_length, si_u, si_nu))));
	LBM lbm(lbm_N, lbm_nu);
	// ###################################################################################### define geometry ######################################################################################
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(180.0f))*float3x3(float3(0, 0, 1), radians(180.0f));
	Mesh* mesh = read_stl(get_exe_path()+"../stl/Cow_t.stl", lbm.size(), lbm.center(), rotation, lbm_length); // https://www.thingiverse.com/thing:182114/files
	mesh->translate(float3(0.0f, 1.0f-mesh->pmin.y+0.1f*lbm_length, 1.0f-mesh->pmin.z)); // move mesh forward a bit and to simulation box bottom, keep in mind 1 cell thick box boundaries
	lbm.voxelize_mesh_on_device(mesh);
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z==0u) lbm.flags[n] = TYPE_S; // solid floor
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u; // initialize y-velocity everywhere except in solid cells
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all other simulation box boundaries are inflow/outflow
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
	lbm.graphics.set_camera_centered(-40.0f, 20.0f, 78.0f, 1.25f);
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 10.0f)) lbm.graphics.write_frame();
		lbm.run(1u, lbm_T);
	}
#else // GRAPHICS && !INTERACTIVE_GRAPHICS
	lbm.run();
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
} /**/



/*void main_setup() { // Space Shuttle; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 4.0f, 0.8f), 1000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_Re = 10000000.0f;
	const float lbm_u = 0.075f;
	const ulong lbm_T = 108000ull;
	LBM lbm(lbm_N, 2u, 4u, 1u, units.nu_from_Re(lbm_Re, (float)lbm_N.x, lbm_u)); // run on 2x4x1 = 8 GPUs
	// ###################################################################################### define geometry ######################################################################################
	const float size = 1.25f*lbm.size().x;
	const float3 center = float3(lbm.center().x, 0.55f*size, lbm.center().z+0.05f*size);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-20.0f))*float3x3(float3(0, 0, 1), radians(270.0f));
	Clock clock;
	lbm.voxelize_stl(get_exe_path()+"../stl/Full_Shuttle.stl", center, rotation, size); // https://www.thingiverse.com/thing:4975964/files
	println(print_time(clock.stop()));
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_FLAG_SURFACE|VIS_Q_CRITERION;
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
	lbm.write_status();
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 30.0f)) {
			lbm.graphics.set_camera_free(float3(-1.435962f*(float)Nx, 0.364331f*(float)Ny, 1.344426f*(float)Nz), -205.0f, 36.0f, 74.0f); // top
			lbm.graphics.write_frame(get_exe_path()+"export/top/");
			lbm.graphics.set_camera_free(float3(-1.021207f*(float)Nx, -0.518006f*(float)Ny, 0.0f*(float)Nz), -137.0f, 0.0f, 74.0f); // bottom
			lbm.graphics.write_frame(get_exe_path()+"export/bottom/");
		}
		lbm.run(1u, lbm_T);
	}
	lbm.write_status();
#else // GRAPHICS && !INTERACTIVE_GRAPHICS
	lbm.run();
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
} /**/



/*void main_setup() { // Starship; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 2.0f), 1000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_Re = 10000000.0f;
	const float lbm_u = 0.05f;
	const ulong lbm_T = 108000ull;
	LBM lbm(lbm_N, 1u, 1u, 1u, units.nu_from_Re(lbm_Re, (float)lbm_N.x, lbm_u));
	// ###################################################################################### define geometry ######################################################################################
	const float size = 1.6f*lbm.size().x;
	const float3 center = float3(lbm.center().x, lbm.center().y+0.05f*size, 0.18f*size);
	lbm.voxelize_stl(get_exe_path()+"../stl/StarShipV2.stl", center, size); // https://www.thingiverse.com/thing:4912729/files
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.z[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_FLAG_SURFACE|VIS_Q_CRITERION;
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
	lbm.write_status();
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 20.0f)) {
			lbm.graphics.set_camera_free(float3(2.116744f*(float)Nx, -0.775261f*(float)Ny, 1.026577f*(float)Nz), -38.0f, 37.0f, 60.0f); // top
			lbm.graphics.write_frame(get_exe_path()+"export/top/");
			lbm.graphics.set_camera_free(float3(0.718942f*(float)Nx, 0.311263f*(float)Ny, -0.498366f*(float)Nz), 32.0f, -40.0f, 104.0f); // bottom
			lbm.graphics.write_frame(get_exe_path()+"export/bottom/");
			lbm.graphics.set_camera_free(float3(1.748119f*(float)Nx, 0.442782f*(float)Ny, 0.087945f*(float)Nz), 24.0f, 2.0f, 92.0f); // side
			lbm.graphics.write_frame(get_exe_path()+"export/side/");
		}
		lbm.run(1u, lbm_T);
	}
	lbm.write_status();
#else // GRAPHICS && !INTERACTIVE_GRAPHICS
	lbm.run();
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
} /**/



/*void main_setup() { // Ahmed body; required extensions in defines.hpp: FP16C, FORCE_FIELD, EQUILIBRIUM_BOUNDARIES, SUBGRID, optionally INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint memory = 10000u; // available VRAM of GPU(s) in MB
	const float lbm_u = 0.05f;
	const float box_scale = 6.0f;
	const float si_u = 60.0f;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	const float si_width=0.389f, si_height=0.288f, si_length=1.044f;
	const float si_A = si_width*si_height+2.0f*0.05f*0.03f;
	const float si_T = 0.25f;
	const float si_Lx = units.x(box_scale*si_width);
	const float si_Ly = units.x(box_scale*si_length);
	const float si_Lz = units.x(0.5f*(box_scale-1.0f)*si_width+si_height);
	const uint3 lbm_N = resolution(float3(si_Lx, si_Ly, si_Lz), memory); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	units.set_m_kg_s((float)lbm_N.y, lbm_u, 1.0f, box_scale*si_length, si_u, si_rho);
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	const float lbm_length = units.x(si_length);
	print_info("Re = "+to_string(to_uint(units.si_Re(si_width, si_u, si_nu))));
	LBM lbm(lbm_N, lbm_nu);
	// ###################################################################################### define geometry ######################################################################################
	Mesh* mesh = read_stl(get_exe_path()+"../stl/ahmed_25deg_m.stl", lbm.size(), lbm.center(), float3x3(float3(0, 0, 1), radians(90.0f)), lbm_length);
	mesh->translate(float3(0.0f, units.x(0.5f*(0.5f*box_scale*si_length-si_width))-mesh->pmin.y, 1.0f-mesh->pmin.z));
	lbm.voxelize_mesh_on_device(mesh, TYPE_S|TYPE_X); // https://github.com/nathanrooy/ahmed-bluff-body-cfd/blob/master/geometry/ahmed_25deg_m.stl converted to binary
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z==0u) lbm.flags[n] = TYPE_S;
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u) lbm.flags[n] = TYPE_E;
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_FIELD;
	lbm.graphics.field_mode = 1;
	lbm.graphics.slice_mode = 1;
	//lbm.graphics.set_camera_centered(20.0f, 30.0f, 10.0f, 1.648722f);
	lbm.run(0u, lbm_T); // initialize simulation
#if defined(FP16S)
	const string path = get_exe_path()+"FP16S/"+to_string(memory)+"MB/";
#elif defined(FP16C)
	const string path = get_exe_path()+"FP16C/"+to_string(memory)+"MB/";
#else // FP32
	const string path = get_exe_path()+"FP32/"+to_string(memory)+"MB/";
#endif // FP32
	//lbm.write_status(path);
	//write_file(path+"Cd.dat", "# t\tCd\n");
	const float3 lbm_com = lbm.object_center_of_mass(TYPE_S|TYPE_X);
	print_info("com = "+to_string(lbm_com.x, 2u)+", "+to_string(lbm_com.y, 2u)+", "+to_string(lbm_com.z, 2u));
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		Clock clock;
		const float3 lbm_force = lbm.object_force(TYPE_S|TYPE_X);
		//const float3 lbm_torque = lbm.object_torque(lbm_com, TYPE_S|TYPE_X);
		//print_info("F="+to_string(lbm_force.x, 2u)+","+to_string(lbm_force.y, 2u)+","+to_string(lbm_force.z, 2u)+", T="+to_string(lbm_torque.x, 2u)+","+to_string(lbm_torque.y, 2u)+","+to_string(lbm_torque.z, 2u)+", t="+to_string(clock.stop(), 3u));
		const float Cd = units.si_F(lbm_force.y)/(0.5f*si_rho*sq(si_u)*si_A); // expect Cd to be too large by a factor 1.3-2.0x; need wall model
		print_info("Cd = "+to_string(Cd, 3u)+", t = "+to_string(clock.stop(), 3u));
		//write_line(path+"Cd.dat", to_string(lbm.get_t())+"\t"+to_string(Cd, 3u)+"\n");
		lbm.run(1u, lbm_T);
	}
	//lbm.write_status(path);
} /**/



/*void main_setup() { // Cessna 172 propeller aircraft; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, MOVING_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 0.8f, 0.25f), 8000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_u = 0.075f;
	const float lbm_width = 0.95f*(float)lbm_N.x;
	const ulong lbm_dt = 4ull; // revoxelize rotor every dt time steps
	const float si_T = 1.0f;
	const float si_width = 11.0f;
	const float si_u = 226.0f/3.6f;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	units.set_m_kg_s(lbm_width, lbm_u, 1.0f, si_width, si_u, si_rho);
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	print_info("Re = "+to_string(to_uint(units.si_Re(si_width, si_u, si_nu))));
	print_info(to_string(si_T, 3u)+" seconds = "+to_string(lbm_T)+" time steps");
	LBM lbm(lbm_N, units.nu(si_nu));
	// ###################################################################################### define geometry ######################################################################################
	Mesh* plane = read_stl(get_exe_path()+"../stl/Cessna-172-Skyhawk-body.stl"); // https://www.thingiverse.com/thing:814319/files
	Mesh* rotor = read_stl(get_exe_path()+"../stl/Cessna-172-Skyhawk-rotor.stl"); // plane and rotor separated with Microsoft 3D Builder
	const float scale = lbm_width/plane->get_bounding_box_size().x; // scale plane and rotor to simulation box size
	plane->scale(scale);
	rotor->scale(scale);
	const float3 offset = lbm.center()-plane->get_bounding_box_center(); // move plane and rotor to simulation box center
	plane->translate(offset);
	rotor->translate(offset);
	plane->set_center(plane->get_center_of_mass()); // set rotation center of mesh to its center of mass
	rotor->set_center(rotor->get_center_of_mass());
	const float lbm_radius=0.5f*rotor->get_max_size(), omega=-lbm_u/lbm_radius, domega=omega*(float)lbm_dt;
	lbm.voxelize_mesh_on_device(plane);
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		lbm.voxelize_mesh_on_device(rotor, TYPE_S, rotor->get_center(), float3(0.0f), float3(0.0f, omega, 0.0f)); // revoxelize mesh on GPU
		lbm.run(lbm_dt, lbm_T); // run dt time steps
		rotor->rotate(float3x3(float3(0.0f, 1.0f, 0.0f), domega)); // rotate mesh
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		if(lbm.graphics.next_frame(lbm_T, 5.0f)) {
			lbm.graphics.set_camera_free(float3(0.192778f*(float)Nx, -0.669183f*(float)Ny, 0.657584f*(float)Nz), -77.0f, 27.0f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/a/");
			lbm.graphics.set_camera_free(float3(0.224926f*(float)Nx, -0.594332f*(float)Ny, -0.277894f*(float)Nz), -65.0f, -14.0f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/b/");
			lbm.graphics.set_camera_free(float3(-0.000000f*(float)Nx, 0.650189f*(float)Ny, 1.461048f*(float)Nz), 90.0f, 40.0f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/c/");
		}
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
	}
} /**/



/*void main_setup() { // Bell 222 helicopter; required extensions in defines.hpp: FP16C, EQUILIBRIUM_BOUNDARIES, MOVING_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 1.2f, 0.3f), 8000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_u = 0.16f;
	const float lbm_length = 0.8f*(float)lbm_N.x;
	const float si_T = 0.34483f; // 2 revolutions of the main rotor
	const ulong lbm_dt = 4ull; // revoxelize rotor every dt time steps
	const float si_length=12.85f, si_d=12.12f, si_rpm=348.0f;
	const float si_u = si_rpm/60.0f*si_d*pif;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	units.set_m_kg_s(lbm_length, lbm_u, 1.0f, si_length, si_u, si_rho);
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	LBM lbm(lbm_N, 1u, 1u, 1u, lbm_nu);
	// ###################################################################################### define geometry ######################################################################################
	Mesh* body = read_stl(get_exe_path()+"../stl/Bell-222-body.stl"); // https://www.thingiverse.com/thing:1625155/files
	Mesh* main = read_stl(get_exe_path()+"../stl/Bell-222-main.stl"); // body and rotors separated with Microsoft 3D Builder
	Mesh* back = read_stl(get_exe_path()+"../stl/Bell-222-back.stl");
	const float scale = lbm_length/body->get_bounding_box_size().y; // scale body and rotors to simulation box size
	body->scale(scale);
	main->scale(scale);
	back->scale(scale);
	const float3 offset = lbm.center()-body->get_bounding_box_center(); // move body and rotors to simulation box center
	body->translate(offset);
	main->translate(offset);
	back->translate(offset);
	body->set_center(body->get_center_of_mass()); // set rotation center of mesh to its center of mass
	main->set_center(main->get_center_of_mass());
	back->set_center(back->get_center_of_mass());
	const float main_radius=0.5f*main->get_max_size(), main_omega=lbm_u/main_radius, main_domega=main_omega*(float)lbm_dt;
	const float back_radius=0.5f*back->get_max_size(), back_omega=-lbm_u/back_radius, back_domega=back_omega*(float)lbm_dt;
	lbm.voxelize_mesh_on_device(body);
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] =  0.2f*lbm_u;
		if(lbm.flags[n]!=TYPE_S) lbm.u.z[n] = -0.1f*lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		lbm.voxelize_mesh_on_device(main, TYPE_S, main->get_center(), float3(0.0f), float3(0.0f, 0.0f, main_omega)); // revoxelize mesh on GPU
		lbm.voxelize_mesh_on_device(back, TYPE_S, back->get_center(), float3(0.0f), float3(back_omega, 0.0f, 0.0f)); // revoxelize mesh on GPU
		lbm.run(lbm_dt, lbm_T); // run dt time steps
		main->rotate(float3x3(float3(0.0f, 0.0f, 1.0f), main_domega)); // rotate mesh
		back->rotate(float3x3(float3(1.0f, 0.0f, 0.0f), back_domega)); // rotate mesh
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		if(lbm.graphics.next_frame(lbm_T, 10.0f)) {
			lbm.graphics.set_camera_free(float3(0.528513f*(float)Nx, 0.102095f*(float)Ny, 1.302283f*(float)Nz), 16.0f, 47.0f, 96.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/a/");
			lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, -0.114244f*(float)Ny, 0.543265f*(float)Nz), 90.0f+degrees((float)lbm.get_t()/(float)lbm_dt*main_domega), 36.0f, 120.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/b/");
			lbm.graphics.set_camera_free(float3(0.557719f*(float)Nx, -0.503388f*(float)Ny, -0.591976f*(float)Nz), -43.0f, -21.0f, 75.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/c/");
			lbm.graphics.set_camera_centered(58.0f, 9.0f, 88.0f, 1.648722f);
			lbm.graphics.write_frame(get_exe_path()+"export/d/");
			lbm.graphics.set_camera_centered(0.0f, 90.0f, 100.0f, 1.100000f);
			lbm.graphics.write_frame(get_exe_path()+"export/e/");
			lbm.graphics.set_camera_free(float3(0.001612f*(float)Nx, 0.523852f*(float)Ny, 0.992613f*(float)Nz), 90.0f, 37.0f, 94.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/f/");
		}
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
	}
} /**/



/*void main_setup() { // Mercedes F1 W14 car; required extensions in defines.hpp: FP16S, EQUILIBRIUM_BOUNDARIES, MOVING_BOUNDARIES, SUBGRID, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 0.5f), 4000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_u = 0.075f;
	const float lbm_length = 0.8f*(float)lbm_N.y;
	const float si_T = 0.25f;
	const float si_u = 100.0f/3.6f;
	const float si_length=5.5f, si_width=2.0f;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	units.set_m_kg_s(lbm_length, lbm_u, 1.0f, si_length, si_u, si_rho);
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	print_info("Re = "+to_string(to_uint(units.si_Re(si_width, si_u, si_nu))));
	LBM lbm(lbm_N, 1u, 1u, 1u, lbm_nu);
	// ###################################################################################### define geometry ######################################################################################
	Mesh* body = read_stl(get_exe_path()+"../stl/mercedesf1-body.stl"); // https://downloadfree3d.com/3d-models/vehicles/sports-car/mercedes-f1-w14/
	Mesh* front_wheels = read_stl(get_exe_path()+"../stl/mercedesf1-front-wheels.stl"); // wheels separated, decals removed and converted to .stl in Microsoft 3D Builder
	Mesh* back_wheels = read_stl(get_exe_path()+"../stl/mercedesf1-back-wheels.stl"); // to avoid instability from too small gaps: remove front wheel fenders and move out right back wheel a bit
	const float scale = lbm_length/body->get_bounding_box_size().y; // scale parts
	body->scale(scale);
	front_wheels->scale(scale);
	back_wheels->scale(scale);
	const float3 offset = float3(lbm.center().x-body->get_bounding_box_center().x, 1.0f-body->pmin.y+0.25f*back_wheels->get_min_size(), 4.0f-back_wheels->pmin.z);
	body->translate(offset);
	front_wheels->translate(offset);
	back_wheels->translate(offset);
	body->set_center(body->get_center_of_mass()); // set rotation center of mesh to its center of mass
	front_wheels->set_center(front_wheels->get_center_of_mass());
	back_wheels->set_center(back_wheels->get_center_of_mass());
	const float lbm_radius=0.5f*back_wheels->get_min_size(), omega=lbm_u/lbm_radius;
	lbm.voxelize_mesh_on_device(body);
	lbm.voxelize_mesh_on_device(front_wheels, TYPE_S, front_wheels->get_center(), float3(0.0f), float3(omega, 0.0f, 0.0f)); // make wheels rotating
	lbm.voxelize_mesh_on_device(back_wheels, TYPE_S, back_wheels->get_center(), float3(0.0f), float3(omega, 0.0f, 0.0f)); // make wheels rotating
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u) lbm.flags[n] = TYPE_E;
		if(z==0u) lbm.flags[n] = TYPE_S;
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 30.0f)) {
			lbm.graphics.set_camera_free(float3(0.779346f*(float)Nx, -0.315650f*(float)Ny, 0.329444f*(float)Nz), -27.0f, 19.0f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/a/");
			lbm.graphics.set_camera_free(float3(0.556877f*(float)Nx, 0.228191f*(float)Ny, 1.159613f*(float)Nz), 19.0f, 53.0f, 100.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/b/");
			lbm.graphics.set_camera_free(float3(0.220650f*(float)Nx, -0.589529f*(float)Ny, 0.085407f*(float)Nz), -72.0f, 16.0f, 86.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/c/");
			const float progress = (float)lbm.get_t()/(float)lbm_T;
			const float A = 75.0f, B = -160.0f;
			lbm.graphics.set_camera_centered(A+progress*(B-A), -5.0f, 100.0f, 1.648721f);
			lbm.graphics.write_frame(get_exe_path()+"export/d/");
		}
		lbm.run(1u, lbm_T);
	}
#else // GRAPHICS && !INTERACTIVE_GRAPHICS
	lbm.run();
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
} /**/



/*void main_setup() { // hydraulic jump; required extensions in defines.hpp: FP16S, VOLUME_FORCE, EQUILIBRIUM_BOUNDARIES, MOVING_BOUNDARIES, SURFACE, SUBGRID, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint memory = 208u; // GPU VRAM in MB
	const float si_T = 100.0f; // simulated time in [s]

	const float3 si_N = float3(0.96f, 3.52f, 0.96f); // box size in [m]
	const float si_p1 = si_N.y*3.0f/20.0f; // socket length in [m]
	const float si_h1 = si_N.z*2.0f/5.0f; // socket height in [m]
	const float si_h2 = si_N.z*3.0f/5.0f; // water height in [m]

	const float si_Q = 0.25f; // inlet volumetric flow rate in [m^3/s]
	const float si_A_inlet = si_N.x*(si_h2-si_h1); // inlet cross-section area in [m^2]
	const float si_A_outlet = si_N.x*si_h1; // outlet cross-section area in [m^2]
	const float si_u_inlet = si_Q/si_A_inlet; // inlet average flow velocity in [m/s]
	const float si_u_outlet = si_Q/si_A_outlet; // outlet average flow velocity in [m/s]

	float const si_nu = 1.0E-6f; // kinematic shear viscosity [m^2/s]
	const float si_rho = 1000.0f; // water density [kg/m^3]
	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
	//const float si_sigma = 73.81E-3f; // water surface tension [kg/s^2] (no need to use surface tension here)

	const uint3 lbm_N = resolution(si_N, memory); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_u_inlet = 0.075f; // velocity in LBM units for pairing lbm_u with si_u --> lbm_u in LBM units will be equivalent si_u in SI units
	units.set_m_kg_s((float)lbm_N.y, lbm_u_inlet, 1.0f, si_N.y, si_u_inlet, si_rho); // calculate 3 independent conversion factors (m, kg, s)

	const float lbm_nu = units.nu(si_nu); // kinematic shear viscosity
	const ulong lbm_T = units.t(si_T); // how many time steps to compute to cover exactly si_T seconds in real time
	const float lbm_f = units.f(si_rho, si_g); // force per volume
	//const float lbm_sigma = units.sigma(si_sigma); // surface tension (not required here)

	const uint lbm_p1 = to_uint(units.x(si_p1));
	const uint lbm_h1 = to_uint(units.x(si_h1));
	const uint lbm_h2 = to_uint(units.x(si_h2));
	const float lbm_u_outlet = units.u(si_u_outlet);

	LBM lbm(lbm_N, 1u, 1u, 1u, lbm_nu, 0.0f, 0.0f, -lbm_f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z<lbm_h2) {
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(0.0005f, z, lbm_h2);
		}
		if(y<lbm_p1&&z<lbm_h1) lbm.flags[n] = TYPE_S;
		if(y<=1u&&x>0u&&x<Nx-1u&&z>=lbm_h1&&z<lbm_h2) {
			lbm.flags[n] = y==0u ? TYPE_S : TYPE_F;
			lbm.u.y[n] = lbm_u_inlet;
		}
		if(y==Ny-1u&&x>0u&&x<Nx-1u&&z>0u) {
			lbm.flags[n] = TYPE_E;
			lbm.u.y[n] = lbm_u_outlet;
		}
		if(x==0u||x==Nx-1u||y==0u||z==0u) lbm.flags[n] = TYPE_S; // sides and bottom non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
	lbm.run();
	//lbm.run(1000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/2u, Ny/4u, Nz/4u)]); wait(); // test for binary identity
} /**/



/*void main_setup() { // dam break; required extensions in defines.hpp: FP16S, VOLUME_FORCE, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(128u, 256u, 256u, 0.005f, 0.0f, 0.0f, -0.0002f, 0.0001f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z<Nz*6u/8u && y<Ny/8u) lbm.flags[n] = TYPE_F;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
	lbm.run();
} /**/



/*void main_setup() { // liquid metal on a speaker; required extensions in defines.hpp: FP16S, VOLUME_FORCE, MOVING_BOUNDARIES, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint L = 128u;
	const float u = 0.09f; // peak velocity of speaker membrane
	const float f = 0.0005f;
	const float frequency = 0.01f; // amplitude = u/(2.0f*pif*frequency);
	LBM lbm(L, L, L*3u/4u, 0.01f, 0.0f, 0.0f, -f, 0.005f);
	// ###################################################################################### define geometry ######################################################################################
	const uint threads = (uint)thread::hardware_concurrency();
	vector<uint> seed(threads);
	for(uint t=0u; t<threads; t++) seed[t] = 42u+t;
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), threads, [&](ulong n, uint t) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z<Nz/3u && x>0u&&x<Nx-1u&&y>0u&&y<Ny-1u&&z>0u&&z<Nz-1u) {
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)(Nz/3u));
			lbm.u.x[n] = random_symmetric(seed[t], 1E-9f);
			lbm.u.y[n] = random_symmetric(seed[t], 1E-9f);
			lbm.u.z[n] = random_symmetric(seed[t], 1E-9f);
			lbm.flags[n] = TYPE_F;
		}
		if(z==0u) lbm.u.z[n] = 1E-16f;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
	lbm.run(0u); // initialize simulation
	while(true) { // main simulation loop
		lbm.u.read_from_device();
		const float uz = u*sinf(2.0f*pif*frequency*(float)lbm.get_t());
		for(uint z=0u; z<1u; z++) {
			for(uint y=1u; y<Ny-1u; y++) {
				for(uint x=1u; x<Nx-1u; x++) {
					const uint n = x+(y+z*Ny)*Nx;
					lbm.u.z[n] = uz;
				}
			}
		}
		lbm.u.write_to_device();
		lbm.run(1u);
	}
} /**/



/*void main_setup() { // breaking waves on beach; required extensions in defines.hpp: FP16S, VOLUME_FORCE, EQUILIBRIUM_BOUNDARIES, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const float f = 0.001f; // make smaller
	const float u = 0.12f; // peak velocity of speaker membrane
	const float frequency = 0.0007f; // amplitude = u/(2.0f*pif*frequency);
	LBM lbm(128u, 640u, 96u, 0.01f, 0.0f, 0.0f, -f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		const uint H = Nz/2u;
		if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		}
		if(plane(x, y, z, float3(lbm.center().x, 128.0f, 0.0f), float3(0.0f, -1.0f, 8.0f))) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
		if(y==0u && x>0u&&x<Nx-1u&&z>0u&&z<Nz-1u) lbm.flags[n] = TYPE_E;
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE | (lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE);
	lbm.run(0u); // initialize simulation
	while(true) { // main simulation loop
		lbm.u.read_from_device();
		const float uy = u*sinf(2.0f*pif*frequency*(float)lbm.get_t());
		const float uz = 0.5f*u*cosf(2.0f*pif*frequency*(float)lbm.get_t());
		for(uint z=1u; z<Nz-1u; z++) {
			for(uint y=0u; y<1u; y++) {
				for(uint x=1u; x<Nx-1u; x++) {
					const uint n = x+(y+z*Ny)*Nx;
					lbm.u.y[n] = uy;
					lbm.u.z[n] = uz;
				}
			}
		}
		lbm.u.write_to_device();
		lbm.run(100u);
	}
} /**/



/*void main_setup() { // river; required extensions in defines.hpp: FP16S, VOLUME_FORCE, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(128u, 384u, 96u, 0.02f, 0.0f, -0.00007f, -0.0005f, 0.01f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		const int R = 20, H = 32;
		if(z==0) lbm.flags[n] = TYPE_S;
		else if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = -0.1f;
		}
		if(cylinder(x, y, z, float3(Nx*2u/3u, Ny*2u/3u, Nz/2u)+0.5f, float3(0u, 0u, Nz), (float)R)) lbm.flags[n] = TYPE_S;
		if(cuboid(x, y, z, float3(Nx/3u, Ny/3u, Nz/2u)+0.5f, float3(2u*R, 2u*R, Nz))) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u) lbm.flags[n] = TYPE_S; // x non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
	lbm.run();
} /**/



/*void main_setup() { // raindrop impact; required extensions in defines.hpp: FP16C, VOLUME_FORCE, EQUILIBRIUM_BOUNDARIES, SURFACE, INTERACTIVE_GRAPHICS or GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(1.0f, 1.0f, 0.85f), 4000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	float lbm_D = (float)lbm_N.x/5.0f;
	const float lbm_u = 0.05f; // impact velocity in LBM units
	const float si_T = 0.003f; // simulated time in [s]
	const float inclination = 20.0f; // impact angle [°], 0 = vertical
	const int select_drop_size = 12;
	//                            0        1        2        3        4        5        6        7        8        9       10       11       12       13 (13 is for validation)
	const float si_Ds[] = { 1.0E-3f, 1.5E-3f, 2.0E-3f, 2.5E-3f, 3.0E-3f, 3.5E-3f, 4.0E-3f, 4.5E-3f, 5.0E-3f, 5.5E-3f, 6.0E-3f, 6.5E-3f, 7.0E-3f, 4.1E-3f };
	const float si_us[] = {   4.50f,   5.80f,   6.80f,   7.55f,   8.10f,   8.45f,   8.80f,   9.05f,   9.20f,   9.30f,   9.40f,   9.45f,   9.55f,   7.21f };
	float const si_nu = 1.0508E-6f; // kinematic shear viscosity [m^2/s] at 20°C and 35g/l salinity
	const float si_rho = 1024.8103f; // fluid density [kg/m^3] at 20°C and 35g/l salinity
	const float si_sigma = 73.81E-3f; // fluid surface tension [kg/s^2] at 20°C and 35g/l salinity
	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
	const float si_D = si_Ds[select_drop_size]; // drop diameter [m] (1-7mm)
	const float si_u = si_us[select_drop_size]; // impact velocity [m/s] (4.50-9.55m/s)
	units.set_m_kg_s(lbm_D, lbm_u, 1.0f, si_D, si_u, si_rho); // calculate 3 independent conversion factors (m, kg, s)
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	const float lbm_f = units.f(si_rho, si_g);
	const float lbm_sigma = units.sigma(si_sigma);
	print_info("D = "+to_string(si_D, 6u));
	print_info("Re = "+to_string(units.si_Re(si_D, si_u, si_nu), 6u));
	print_info("We = "+to_string(units.si_We(si_D, si_u, si_rho, si_sigma), 6u));
	print_info("Fr = "+to_string(units.si_Fr(si_D, si_u, si_g), 6u));
	print_info("Ca = "+to_string(units.si_Ca(si_u, si_rho, si_nu, si_sigma), 6u));
	print_info("Bo = "+to_string(units.si_Bo(si_D, si_rho, si_g, si_sigma), 6u));
	print_info(to_string(to_uint(1000.0f*si_T))+" ms = "+to_string(units.t(si_T))+" LBM time steps");
	const float lbm_H = 0.4f*(float)lbm_N.x;
	const float lbm_R = 0.5f*lbm_D; // drop radius
	LBM lbm(lbm_N, 1u, 1u, 1u, lbm_nu, 0.0f, 0.0f, -lbm_f, lbm_sigma); // calculate values for remaining parameters in simulation units
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(sphere(x, y, z, float3(0.5f*(float)Nx, 0.5f*(float)Ny-2.0f*lbm_R*tan(inclination*pif/180.0f), lbm_H+lbm_R+2.5f)+0.5f, lbm_R+2.0f)) {
			const float b = sphere_plic(x, y, z, float3(0.5f*(float)Nx, 0.5f*(float)Ny-2.0f*lbm_R*tan(inclination*pif/180.0f)+0.5f, lbm_H+lbm_R+2.5f), lbm_R);
			if(b!=-1.0f) {
				lbm.u.y[n] =  sinf(inclination*pif/180.0f)*lbm_u;
				lbm.u.z[n] = -cosf(inclination*pif/180.0f)*lbm_u;
				if(b==1.0f) {
					lbm.flags[n] = TYPE_F;
					lbm.phi[n] = 1.0f;
				} else {
					lbm.flags[n] = TYPE_I;
					lbm.phi[n] = b; // initialize cell fill level phi directly instead of just flags, this way the raindrop sphere is smooth already at initialization
				}
			}
		}
		if(z==0) lbm.flags[n] = TYPE_S;
		else if(z==to_uint(lbm_H)) {
			lbm.flags[n] = TYPE_I;
			lbm.phi[n] = 0.5f; // not strictly necessary, but should be clearer (phi is automatically initialized to 0.5f for TYPE_I if not initialized)
		} else if((float)z<lbm_H) lbm.flags[n] = TYPE_F;
		else if((x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u)&&(float)z>lbm_H+0.5f*lbm_R) { // make drops that hit the simulation box ceiling disappear
			lbm.rho[n] = 0.5f;
			lbm.flags[n] = TYPE_E;
		}
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS) && !defined(INTERACTIVE_GRAPHICS_ASCII)
	lbm.run(0u, lbm_T); // initialize simulation
	while(lbm.get_t()<=lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 20.0f)) { // generate video
			lbm.graphics.set_camera_centered(-30.0f, 20.0f, 100.0f, 1.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/n/");
			lbm.graphics.set_camera_centered(10.0f, 40.0f, 100.0f, 1.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/p/");
			lbm.graphics.set_camera_centered(0.0f, 0.0f, 45.0f, 1.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/o/");
			lbm.graphics.set_camera_centered(0.0f, 90.0f, 45.0f, 1.0f);
			lbm.graphics.write_frame(get_exe_path()+"export/t/");
		}
		lbm.run(1u, lbm_T);
	}
	//lbm.run(lbm_T); // only generate one image
	//lbm.graphics.set_camera_centered(-30.0f, 20.0f, 100.0f, 1.0f);
	//lbm.graphics.write_frame();
#else // GRAPHICS && !INTERACTIVE_GRAPHICS
	lbm.run();
#endif // GRAPHICS && !INTERACTIVE_GRAPHICS
} /**/



/*void main_setup() { // bursting bubble; required extensions in defines.hpp: FP16C, VOLUME_FORCE, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(4.0f, 4.0f, 3.0f), 1000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float lbm_d = 0.25f*(float)lbm_N.x; // bubble diameter in LBM units
	const float lbm_sigma = 0.0003f; // surface tension coefficient in LBM units
	const float si_nu = 1E-6f; // kinematic shear viscosity (water) [m^2/s]
	const float si_rho = 1E3f; // density (water) [kg/m^3]
	const float si_sigma = 0.072f; // surface tension (water) [kg/s^2]
	const float si_d = 4E-3f; // bubble diameter [m]
	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
	const float si_f = units.si_f_from_si_g(si_g, si_rho);
	const float lbm_rho = 1.0f;
	const float m = si_d/lbm_d; // length si_x = x*[m]
	const float kg = si_rho/lbm_rho*cb(m); // density si_rho = rho*[kg/m^3]
	const float s = sqrt(lbm_sigma/si_sigma*kg); // velocity si_sigma = sigma*[kg/s^2]
	units.set_m_kg_s(m, kg, s); // do unit conversion manually via d, rho and sigma
	const uint lbm_H = to_uint(2.0f*lbm_d);
	LBM lbm(lbm_N, units.nu(si_nu), 0.0f, 0.0f, -units.f(si_f), lbm_sigma);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(z<lbm_H) lbm.flags[n] = TYPE_F;
		const float r = 0.5f*lbm_d;
		if(sphere(x, y, z, float3(lbm.center().x, lbm.center().y, (float)lbm_H-0.5f*lbm_d), r+1.0f)) { // bubble
			const float b = clamp(sphere_plic(x, y, z, float3(lbm.center().x, lbm.center().y, (float)lbm_H-0.5f*lbm_d), r), 0.0f, 1.0f);
			if(b==1.0f) {
				lbm.flags[n] = TYPE_G;
				lbm.phi[n] = 0.0f;
			} else {
				lbm.flags[n] = TYPE_I;
				lbm.phi[n] = (1.0f-b); // initialize cell fill level phi directly instead of just flags, this way the bubble sphere is smooth already at initialization
			}
		}
		if(z==0) lbm.flags[n] = TYPE_S;
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
	lbm.run();
} /**/



/*void main_setup() { // cube with changing gravity; required extensions in defines.hpp: FP16S, VOLUME_FORCE, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(96u, 96u, 96u, 0.02f, 0.0f, 0.0f, -0.001f, 0.001f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(x<Nx*2u/3u&&y<Ny*2u/3u) lbm.flags[n] = TYPE_F;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S;
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
	lbm.run(0u); // initialize simulation
	while(true) { // main simulation loop
		lbm.set_f(0.0f, 0.0f, -0.001f);
		lbm.run(2500u);
		lbm.set_f(0.0f, +0.001f, 0.0f);
		lbm.run(2500u);
		lbm.set_f(0.0f, 0.0f, +0.001f);
		lbm.run(2500u);
		lbm.set_f(0.0f, -0.001f, 0.0f);
		lbm.run(2000u);
		lbm.set_f(0.0f, 0.0f, 0.0f);
		lbm.run(3000u);
	}
} /**/



/*void main_setup() { // periodic faucet mass conservation test; required extensions in defines.hpp: FP16S, VOLUME_FORCE, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(96u, 192u, 128u, 0.02f, 0.0f, 0.0f, -0.00025f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(y>Ny*5u/6u) lbm.flags[n] = TYPE_F;
		const uint D=max(Nx, Nz), R=D/6;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S; // x and y non periodic
		if((z==0u||z==Nz-1u) && sq(x-Nx/2)+sq(y-Nx/2)>sq(R)) lbm.flags[n] = TYPE_S; // z non periodic
		if(y<=Nx/2u+2u*R && torus_x(x, y, z, float3(Nx/2u, Nx/2u+R, Nz)+0.5f, (float)R, (float)R)) lbm.flags[n] = TYPE_S;
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_PHI_RASTERIZE;
	lbm.run();
} /**/



/*void main_setup() { // two colliding droplets in force field; required extensions in defines.hpp: FP16S, VOLUME_FORCE, FORCE_FIELD, SURFACE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(256u, 256u, 128u, 0.014f, 0.0f, 0.0f, 0.0f, 0.0001f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(sphere(x, y, z, lbm.center()-float3(0u, 10u, 0u), 32.0f)) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = 0.025f;
		}
		if(sphere(x, y, z, lbm.center()+float3(30u, 40u, 0u), 12.0f)) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = -0.2f;
		}
		lbm.F.x[n] = -0.001f*lbm.relative_position(n).x;
		lbm.F.y[n] = -0.001f*lbm.relative_position(n).y;
		lbm.F.z[n] = -0.0005f*lbm.relative_position(n).z;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE;
	lbm.run();
} /**/



/*void main_setup() { // Rayleigh-Benard convection; required extensions in defines.hpp: FP16S, VOLUME_FORCE, TEMPERATURE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(256u, 256u, 64u, 1u, 1u, 1u, 0.02f, 0.0f, 0.0f, -0.0005f, 0.0f, 1.0f, 1.0f);
	// ###################################################################################### define geometry ######################################################################################
	const uint threads = (uint)thread::hardware_concurrency();
	vector<uint> seed(threads);
	for(uint t=0u; t<threads; t++) seed[t] = 42u+t;
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), threads, [&](ulong n, uint t) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		lbm.u.x[n] = random_symmetric(seed[t], 0.015f); // initialize velocity with random noise
		lbm.u.y[n] = random_symmetric(seed[t], 0.015f);
		lbm.u.z[n] = random_symmetric(seed[t], 0.015f);
		lbm.rho[n] = units.rho_hydrostatic(0.0005f, (float)z, 0.5f*(float)Nz); // initialize density with hydrostatic pressure
		if(z==1u) {
			lbm.T[n] = 1.75f;
			lbm.flags[n] = TYPE_T;
		} else if(z==Nz-2u) {
			lbm.T[n] = 0.25f;
			lbm.flags[n] = TYPE_T;
		}
		if(z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // leave lateral simulation box walls periodic by not closing them with TYPE_S
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_STREAMLINES;
	lbm.run();
} /**/



/*void main_setup() { // thermal convection; required extensions in defines.hpp: FP16S, VOLUME_FORCE, TEMPERATURE, INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	LBM lbm(32u, 196u, 60u, 1u, 1u, 1u, 0.02f, 0.0f, 0.0f, -0.0005f, 0.0f, 1.0f, 1.0f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(y==1) {
			lbm.T[n] = 1.8f;
			lbm.flags[n] = TYPE_T;
		} else if(y==Ny-2) {
			lbm.T[n] = 0.3f;
			lbm.flags[n] = TYPE_T;
		}
		lbm.rho[n] = units.rho_hydrostatic(0.0005f, (float)z, 0.5f*(float)Nz); // initialize density with hydrostatic pressure
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_STREAMLINES;
	lbm.run();
	//lbm.run(1000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/2u, Ny/2u, Nz/2u)]); wait(); // test for binary identity
} /**/