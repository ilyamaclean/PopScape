#include <Rcpp.h>
#include <cmath>
#include <unordered_map>
#include <queue>
using namespace Rcpp;
// =========================================================================== #
// ~~~~~~~~~~~ Calculate habitat suitability-weighted distances between ~~~~~~~~~~ #
// ~~~~~~~~~~~~~~ habitat patches, averaged across pixels ~~~~~~~~~~~~~~~~~~~~ #
// =========================================================================== #
// [[Rcpp::export]]
IntegerMatrix renumberPatchIds(IntegerMatrix patchIds) {
    int nrow = patchIds.nrow();
    int ncol = patchIds.ncol();
    // Create an unordered_map to store original patch IDs and their new IDs
    std::unordered_map<int, int> idMap;
    int newId = 1;
    // First pass to assign new sequential IDs to each unique patch ID
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            int patchId = patchIds(i, j);
            if (patchId > 0 && idMap.find(patchId) == idMap.end()) {
                // Assign a new sequential ID for each unique patch ID
                idMap[patchId] = newId++;
            }
        }
    }
    // Second pass to replace original IDs with new sequential IDs
    IntegerMatrix result(nrow, ncol);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            int patchId = patchIds(i, j);
            if (patchId > 0) {
                result(i, j) = idMap[patchId];  // Replace with new ID
            }
            else {
                result(i, j) = 0;  // Leave non-positive IDs unchanged
            }
        }
    }
    return result;
}
// Helper function to calculate Euclidean distance between two pixels
double euclidean_distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}
// Calculate max value of IntegerMatrix
int max_with_na(IntegerMatrix mat) {
    int max_val = INT_MIN; // Smallest possible integer value
    int nrow = mat.nrow();
    int ncol = mat.ncol();
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            if (!IntegerVector::is_na(mat(i, j))) { // Check if not NA
                if (mat(i, j) > max_val) {
                    max_val = mat(i, j); // Update max
                }
            }
        }
    }
    return max_val;
}
// [[Rcpp::export]]
NumericMatrix patchdistcpp(IntegerMatrix patches, NumericMatrix popden) {
    int nrow = patches.nrow();
    int ncol = patches.ncol();
    // Map to store coordinates and habitat quality for each patch's pixels
    std::unordered_map<int, std::vector<std::pair<std::pair<int, int>, double>>> patch_pixels;
    // Gather all pixel locations and qualities by patch ID
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            int patch_id = patches(i, j);
            if (patch_id != NA_INTEGER) { // Skip non-habitat cells
                double den = popden(i, j);
                patch_pixels[patch_id].emplace_back(std::make_pair(i, j), den);
            }
        }
    }
    // Prepare output matrix
    int patch_count = patch_pixels.size();
    NumericMatrix distance_matrix(patch_count, patch_count);
    std::vector<int> patch_ids;
    // Store patch IDs for matrix indexing
    for (const auto& entry : patch_pixels) {
        patch_ids.push_back(entry.first);
    }
    // Calculate weighted mean distances between each pair of patches
    for (int i = 0; i < patch_count; ++i) {
        for (int j = i + 1; j < patch_count; ++j) { // Only compute upper triangle (symmetric)
            int patch_id1 = patch_ids[i];
            int patch_id2 = patch_ids[j];
            double total_weighted_distance = 0.0;
            double total_weight = 0.0;
            // Iterate over all pixels in patch_id1 and patch_id2
            for (const auto& pixel1 : patch_pixels[patch_id1]) {
                int x1 = pixel1.first.first;
                int y1 = pixel1.first.second;
                double quality1 = pixel1.second;
                for (const auto& pixel2 : patch_pixels[patch_id2]) {
                    int x2 = pixel2.first.first;
                    int y2 = pixel2.first.second;
                    double quality2 = pixel2.second;
                    double distance = euclidean_distance(x1, y1, x2, y2);
                    double weight = quality1 * quality2;
                    total_weighted_distance += distance * weight;
                    total_weight += weight;
                }
            }
            // Compute mean weighted distance between patches
            double mean_weighted_distance = total_weight > 0 ? total_weighted_distance / total_weight : NA_REAL;
            distance_matrix(i, j) = mean_weighted_distance;
            distance_matrix(j, i) = mean_weighted_distance; // Symmetric matrix
        }
    }
    return distance_matrix;
}
// Calculate carrying-capacity weighted connectivity
// [[Rcpp::export]]
NumericVector calculateConnectivity(NumericVector Ai, NumericVector Oj, 
    NumericMatrix dij, double alpha)
{
    int n = Ai.size(); // Number of patches
    NumericVector Si(n);
    for (int i = 0; i < n; ++i) {
        double sm = 0.0;
        for (int j = 0; j < n; ++j) {
            double val = 0.0;
            if (Oj[j] > 0.0 && i != j) {
                double Dij = exp(-alpha * dij(i, j));
                val = Dij * Ai[j];
            }
            sm = sm + val;
        }
        Si[i] = sm;
    }
    return Si;
}
// =========================================================================== #
// ~~~~~~~~~~~~~ Calculate average patch size and quality ~~~~~~~~~~~~~~~~~~~~ #
// =========================================================================== #
// [[Rcpp::export]]
NumericVector calculatePatchSizes(IntegerMatrix patches) {
    int nrow = patches.nrow();
    int ncol = patches.ncol();
    int n = max_with_na(patches);
    IntegerVector patchSizes(n);
    // Loop over the matrix to count cells for each patch ID
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            int patch_id = patches(i, j);
            // Check if patches(i, j) is NA
            if (IntegerVector::is_na(patch_id)) {
                patch_id = 0; // Set patch_id to 0 if NA
            }
            if (patch_id > 0) { // Ignore non-habitat cells (ID 0)
                patchSizes[patch_id - 1]++;
            }
        }
    }
    NumericVector patchSv = NumericVector(patchSizes);
    return patchSv;
}
// Calculate average patch quality
// [[Rcpp::export]]
NumericVector calculateAveragePatchQuality(IntegerMatrix patches, NumericMatrix qualities) {
    // Check that both matrices have the same number of cells
    if (patches.nrow() != qualities.nrow() || patches.ncol() != qualities.ncol()) {
        stop("The patch ID and quality matrices must have the same dimensions.");
    }
    int nrow = patches.nrow();
    int ncol = patches.ncol();
    int n = max_with_na(patches);
    NumericVector patchQuals(n);
    IntegerVector patchSizes(n);
    // Loop over the matrix to count cells for each patch ID
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            int patch_id = patches(i, j);
            // Check if patches(i, j) is NA
            if (IntegerVector::is_na(patch_id)) {
                patch_id = 0; // Set patch_id to 0 if NA
            }
            if (patch_id > 0) { // Ignore non-habitat cells (ID 0)
                patchSizes[patch_id - 1]++;
                double val = qualities(i, j);
                if (NumericVector::is_na(val)) {
                    val = 0; // Set quality to 0 if NA
                }
                patchQuals[patch_id - 1] = patchQuals[patch_id - 1] + val;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        patchQuals[i] = patchQuals[i] / patchSizes[i];
    }
    return patchQuals;
}
// =========================================================================== #
// ~~~~~~~~~~~~~~~~~~~~~~~~~~ Run metapopulation model ~~~~~~~~~~~~~~~~~~~~~~~ #
// =========================================================================== #
// ~~~ Calculate extinction probability
// [[Rcpp::export]]
NumericVector ExtProb(NumericVector Ai, NumericVector pden, double mu, double x)
{
    int n = Ai.size();
    NumericVector Ei(n);
    for (int i = 0; i < n; ++i) {
        Ei[i] = mu / pow(Ai[i] * pden[i], x);
        if (Ei[i] > 1.0) Ei[i] = 1.0;
    }
    return Ei;
}
// ~~~ Calculate colonisation probability
// [[Rcpp::export]]
NumericVector ColProb(NumericVector Ai, NumericVector Oj, NumericMatrix dij, 
    double alpha, double gamma)
{
    int n = Ai.size();
    NumericVector Si = calculateConnectivity(Ai, Oj, dij, alpha);
    NumericVector Ci(n);
    for (int i = 0; i < n; ++i) {
        Ci[i] = (Si[i] * Si[i]) / (Si[i] * Si[i] + gamma);
    }
    return Ci;
}
// ~~~ Run simple incident function model
// patch = IntegerMatrix of patch IDs derived from terra::patches
// popden = NumericMatrix of population density derived from a SpatRaster
// res = grid cell resolution in metres
// dij = matrix of patch distances
// mu = extinction probability of a patch of unit size 
// x = scaling of extinction risk with patch area
// alpha = inverse of average dispersal distance in kilometres
// gamma = parameter scaling colonisation to connectivity
// timesteps = number of time steps over which to run model
// [[Rcpp::export]]
NumericVector runmetapopmodel(IntegerMatrix patch, NumericMatrix popden, double res,
    NumericMatrix dij, double mu, double x, double alpha, double gamma, int timesteps)
{
    // Calculate patch sizes
    NumericVector Ai = calculatePatchSizes(patch);
    int n = Ai.size();
    // Return patch size in hectares
    for (int i = 0; i < n; ++i) Ai[i] = (Ai[i] * res * res) / (100.0 * 100.0);
    // Calculate average patch quality
    NumericVector den = calculateAveragePatchQuality(patch, popden);
    // Seed initial patches based on population density
    NumericVector Oj(n);
    for (int i = 0; i < n; ++i) {
        if (den[i] > 0.0) Oj[i] = 1;
    }
    // Calculate extinction probability
    NumericVector Ei = ExtProb(Ai, den, mu, x);
    // Initialize variables
    NumericVector Ci(n);
    NumericVector EiR(n);
    for (int t = 0; t < timesteps; ++t) {
        // Calculate colonisation probability
        Ci = ColProb(Ai, Oj, dij, alpha, gamma);
        for (int i = 0; i < n; ++i) EiR[i] = Ei[i] - Ci[i] * Ei[i];
        // Simulate metapopulation processes
        for (int i = 0; i < n; ++i) {
            // Recalculate extinction probability with rescue effect
            double temp = Ci[i];
            if (Oj[i] > 0.0) {
                temp = 1.0 - EiR[i];
            }
            NumericVector v = rbinom(1, 1, temp);
            Oj[i] = v[0];
        }
    }
    return Oj;
}
// Assign occupancies to patch IDs
// [[Rcpp::export]]
NumericMatrix assignPatchValues(NumericMatrix patchIds, NumericVector values) {
    int nrow = patchIds.nrow();
    int ncol = patchIds.ncol();
    NumericMatrix result(nrow, ncol); // Initialize the result matrix
    // Iterate through each cell in the patch ID matrix
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            int patchId = static_cast<int>(patchIds(i, j));
            // If patchId is valid and within range, assign the corresponding value
            if (patchId > 0 && patchId <= values.size()) {
                result(i, j) = values[patchId - 1]; // Subtract 1 for 0-based index
            }
            else {
                result(i, j) = NA_REAL; // Assign NA if patchId is out of range
            }
        }
    }
    return result;
}
// ======================================================================== //
// ~~~ Runs poulation model across landscape patches for one time step ~~~~ //
// ======================================================================== //
// Nt - vector of initial populations
// K - NumericVector of carrying capacities
// birthrate - NumericVector of average number of offpsring
// survival - NumericVector of average survival rate
// dij - NumericMatrix of inter-patch distances
// dispdist - average disperal distance
// fracdisp - fraction of population emmigrating
// emmsurv - fraction of emmigrating individuals surviving
// Areas - NumericVector of patch areas
// [[Rcpp::export]]
IntegerVector PopSim_one(IntegerVector Nt, IntegerVector K,
    NumericVector birthrate, NumericVector survival, NumericMatrix dij,
    double dispdist, double fracdisp, double emmsurv, NumericVector Areas)
{
    // Sumulate within patch population change
    int pops = Nt.size();
    IntegerVector N(pops);
    for (int i = 0; i < pops; ++i) {
        if (Nt[i] > 0) {
            double lambda = birthrate[i] * survival[i];
            double n = Nt[i];
            NumericVector ind = rpois(n, lambda);
            int tot = 0;
            for (int j = 0; j < n; ++j) tot += ind[j];
            double mu = (1.0 - static_cast<double>(Nt[i]) / K[i]);
            double val = (tot - Nt[i]) * mu + Nt[i];
            N[i] = static_cast<int>(val);
        }
    }
    // Simulate immigration and emmigration
    IntegerVector Emm(pops);
    // ~~ Emmigration
    for (int i = 0; i < pops; ++i) {
        Emm[i] = static_cast<int> (N[i] * fracdisp);
        N[i] += Emm[i];
    }
    // ~~ Immigration
    double alpha = 1 / dispdist;
    for (int i = 0; i < pops; ++i) {
        // cross-sectional area in km 
        double csa = sqrt(Areas[i] / 100);
        double toti = 0.0;
        for (int j = 0; j < pops; ++j) {
            if (i != j) {
                // dispersal liklihood of one individual
                double dl = exp(-alpha * dij(i, j));
                // interception fraction
                double intf = csa / (2 * M_PI * dij(i, j));
                if (intf > 1.0) intf = 1.0;
                // cumulative total number of immigrants from each patch
                toti += dl * intf * Emm[j] * emmsurv;
            }
        }
        N[i] = N[i] - Emm[i] + static_cast<int>(toti);
        if (N[i] > K[i]) N[i] = K[i];
        if (N[i] < 0) N[i] = 0;
    }
    return N;
}
// Nt - vector of initial populations in each patch
// K - IntegerMatrix of carrying capacities in each patch and time-step
// birthrate - NumericMatrix of average number of offpsring in each patch and time-step
// survival - NumericMatrix of average survival rate in each patch and time-step
// dij - NumericMatrix of inter-patch distances
// dispdist - average disperal distance
// fracdisp - fraction of population emmigrating
// emmsurv - fraction of emmigrating individuals surviving
// Areas - NumericVector of patch areas
// [[Rcpp::export]]
IntegerMatrix PopSim(IntegerVector Nt, IntegerMatrix K,
    NumericMatrix birthrate, NumericMatrix survival, NumericMatrix dij,
    double dispdist, double fracdisp, double emmsurv, NumericVector Areas)
{
    int patchn = K.nrow();
    int timesteps = K.ncol();
    IntegerMatrix nout(patchn, timesteps);
    for (int t = 0; t < timesteps; ++t) {
        IntegerVector Kv = K.column(t);
        NumericVector birthratev = birthrate.column(t);
        NumericVector survivalv = survival.column(t);
        Nt = PopSim_one(Nt, Kv, birthratev, survivalv,
            dij, dispdist, fracdisp, emmsurv, Areas);
        for (int i = 0; i < patchn; ++i) nout(i, t) = Nt[i];
    }
    return nout;
}
// Turn a 3D array derived from a StackedRaster into a NumericMatrix by patch
// [[Rcpp::export]]
NumericMatrix arraytomat(NumericVector a, IntegerMatrix patchm)
{
    // Get the dimensions of the 3D array
    IntegerVector dims = a.attr("dim");
    int nrows = dims[0];
    int ncols = dims[1];
    int ntime = dims[2];
    // Find the maximum value in the cleaned vector
    IntegerVector pv = na_omit(as<IntegerVector>(patchm));
    int maxid = max(pv);
    // Initialize the result matrix
    NumericMatrix abundance(maxid, ntime); // Patches are 0-indexed
    IntegerMatrix cell_count(maxid, ntime); // Count of valid cells
    // Iterate through the 3D array
    for (int t = 0; t < ntime; ++t) {
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < ncols; ++j) {
                int patch_id = patchm(i, j);
                if (!IntegerMatrix::is_na(patch_id)) { // Ignore non-patch cells (e.g., NA as -1)
                    abundance(patch_id - 1, t) += a[i + nrows * (j + ncols * t)];
                    cell_count(patch_id - 1, t) += 1;
                }
            }
        }
    }
    // Compute the average
    NumericMatrix result(maxid, ntime);
    for (int patch = 0; patch < maxid; ++patch) {
        for (int t = 0; t < ntime; ++t) {
            if (cell_count(patch, t) > 0) {
                result(patch, t) = abundance(patch, t) / cell_count(patch, t);
            }
            else {
                result(patch, t) = NA_REAL; // No valid cells
            }
        }
    }
    return result;
}
// [[Rcpp::export]]
NumericMatrix apply3D(NumericVector a) {
    // Get the dimensions of the 3D array
    IntegerVector dims = a.attr("dim");
    int nrows = dims[0];
    int ncols = dims[1];
    int ntime = dims[2];
    // Initialize the result matrix
    NumericMatrix result(nrows, ncols);
    // Iterate through each pixel
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            double sum = 0.0;
            int count = 0;
            // Iterate through the time dimension
            for (int t = 0; t < ntime; ++t) {
                double value = a[i + nrows * (j + ncols * t)];
                if (!NumericVector::is_na(value)) { // Skip NA values
                    sum += value;
                    count++;
                }
            }
            // Compute the mean, assign NA if no valid values
            if (count > 0) {
                result(i, j) = sum / count;
            }
            else {
                result(i, j) = NA_REAL;
            }
        }
    }
    return result;
}
IntegerVector aperm3D(IntegerVector a, int rows, int cols, int tsteps) {
    // Initialize a new vector to store the permuted array
    IntegerVector permuted(a.size());
    // Permute the dimensions
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < tsteps; ++k) {
                // Calculate the index in the original array
                int index_original = k + tsteps * (j + cols * i);
                // Calculate the index in the permuted array
                int index_permuted = i + rows * (j + cols * k);
                // Copy the element from the original array to the permuted array
                permuted[index_permuted] = a[index_original];
            }
        }
    }
    // Set the dimensions of the permuted array
    permuted.attr("dim") = IntegerVector::create(rows, cols, tsteps);
    return permuted;
}
// [[Rcpp::export]]
IntegerVector mattoarray(IntegerMatrix m, IntegerMatrix patchm)
{
    // get rows, colums and timesteps
    int rows = patchm.nrow();
    int cols = patchm.ncol();
    int timesteps = m.ncol();
    IntegerVector result(rows * cols * timesteps);
    for (int t = 0; t < timesteps; ++t) {
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                int patch_id = patchm(j, i);
                if (!IntegerMatrix::is_na(patch_id)) { // Ignore non-patch cells (e.g., NA as -1)
                    result[j + rows * (i + cols * t)] = m(patch_id - 1, t);
                }
            }
        }
    }
    // Convert to 3D array
    result.attr("dim") = IntegerVector::create(rows, cols, timesteps);
    //result = aperm3D(result, rows, cols, timesteps);
    return result;
}
// [[Rcpp::export]]
NumericVector apply3Dv(NumericVector a) {
    // Get the dimensions of the 3D array
    IntegerVector dims = a.attr("dim");
    int nrows = dims[0];
    int ncols = dims[1];
    int ntime = dims[2];
    // Initialize the result vector
    NumericVector result(ntime);
    // Iterate through each pixel
    for (int t = 0; t < ntime; ++t) {
        double sum = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < ncols; ++j) {
                double value = a[i + nrows * (j + ncols * t)];
                if (!NumericVector::is_na(value)) { // Skip NA values
                    sum += value;
                }
            }
        }
        result[t] = sum;
    }
    return result;
}