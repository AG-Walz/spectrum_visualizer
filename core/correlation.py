"""R² correlation metrics for spectrum similarity quantification.

Calculates R² (coefficient of determination) between spectra using linear regression.
Provides spectrum_correlation() and new_spectrum_correlation() with Pearson/Cosine methods.
Used to compare experimental vs predicted or experimental vs experimental spectra.
"""
import sys; from pathlib import Path as P; sys.path.insert(0, str(P(__file__).parent.parent / "src"))
from shared_imports import np, pd, Tuple, defaultdict

try:
    from sklearn.linear_model import LinearRegression
except ImportError:
    raise ImportError("scikit-learn is required for R² calculations but not available")


def spectrum_correlation(
    matching_ions_target1: tuple, matching_ions_target2: tuple
) -> tuple:
    """
    Fit linear regression model onto overlapping matching ions between two spectra.
    Process two tuple of (x, y, ion_label) and compute the goodness of fit (R2).
    :param matching_ions_target1: tuple of x, y, ion_label
    :param matching_ions_target2: tuple of x, y, ion_label
    :return: tuple of Goodness of fit (R2), predicted y values, training x and y values
    """
    # Strip ion mass and get overlapping ion labels
    overlapping_ions = set([ion[2] for ion in matching_ions_target1]) & set(
        [ion[2] for ion in matching_ions_target2]
    )
    if len(overlapping_ions) == 0:
        return 0, None, None, None

    # Create dictionaries for easy lookup by ion label
    target1_dict = {ion[2]: ion[1] for ion in matching_ions_target1 if ion[2] in overlapping_ions}
    target2_dict = {ion[2]: ion[1] * -1 for ion in matching_ions_target2 if ion[2] in overlapping_ions}
    
    # Sort ion labels to ensure consistent ordering
    sorted_ions = sorted(overlapping_ions)
    
    # Build arrays in the same order
    target1_intensities = np.array([target1_dict[ion] for ion in sorted_ions])
    target2_intensities = np.array([target2_dict[ion] for ion in sorted_ions])
    
    # Ensure we have enough data points
    if len(target1_intensities) < 2:
        return 0, None, None, None
    
    # Fit a linear regression model
    try:
        reg = LinearRegression().fit(
            target1_intensities.reshape(-1, 1), target2_intensities
        )

        # predict y values using the model
        y_pred = reg.predict(target1_intensities.reshape(-1, 1))

        # calculate R-squared (coefficient of determination) as a measure of goodness of fit
        r2 = reg.score(target1_intensities.reshape(-1, 1), target2_intensities)

        return r2, y_pred, target1_intensities, target2_intensities
    except Exception as e:
        # Log the error and return 0
        import logging
        logger = logging.getLogger(__name__)
        logger.error(f"LinearRegression failed: {e}")
        return 0, None, None, None


def new_spectrum_correlation(
    all_peaks_target1: pd.DataFrame, all_peaks_target2: pd.DataFrame, bins: int = 1000
) -> tuple:
    """
    Fit linear regression model onto binned intensities between two spectra.
    Process two DataFrames and compute the goodness of fit (R2).
    :param all_peaks_target1: DataFrame with columns including 'mz' and 'intensity'
    :param all_peaks_target2: DataFrame with columns including 'mz' and 'intensity'
    :param bins: Number of bins for intensity binning
    :return: tuple of Goodness of fit (R2), predicted y values, training x and y values
    """

    def bin_intensity(df: pd.DataFrame, bins: int, min_mz: float, max_mz: float) -> dict:
        """
        Iterate over m/z according to bin size, catch the peak values for each bin and take
        max -> https://pubs.acs.org/doi/10.1021/ac102272r
        or sum up -> https://link.springer.com/content/pdf/10.1007/978-1-4939-6740-7_7.pdf
        """
        bin_dict = defaultdict(list)
        
        # Use 'mz' column name to match SpecViz DataFrame structure
        mz_col = 'mz' if 'mz' in df.columns else 'Experimental_mass'
        intensity_col = 'intensity' if 'intensity' in df.columns else 'Relative_Intensity'
        
        bin_idx = np.digitize(
            df[mz_col].values, np.linspace(min_mz, max_mz, bins)
        )
        # Create a dictionary, where the keys are the bins and the values are an array of binned mz values
        for b, intensity in zip(bin_idx, df[intensity_col].values):
            bin_dict[b].append(intensity)

        sum_bins = {key: sum(value) for key, value in dict(bin_dict).items()}

        return sum_bins

    # Determine column names based on DataFrame structure
    mz_col1 = 'mz' if 'mz' in all_peaks_target1.columns else 'Experimental_mass'
    mz_col2 = 'mz' if 'mz' in all_peaks_target2.columns else 'Experimental_mass'
    
    # Get bin intensities
    min_mz = round(
        min(
            all_peaks_target1[mz_col1].min(),
            all_peaks_target2[mz_col2].min(),
        ),
        0,
    )
    max_mz = round(
        max(
            all_peaks_target1[mz_col1].max(),
            all_peaks_target2[mz_col2].max(),
        ),
        0,
    )
    bins_target1 = bin_intensity(all_peaks_target1, bins, min_mz, max_mz)
    bins_target2 = bin_intensity(all_peaks_target2, bins, min_mz, max_mz)

    # Get intersection of bins
    bin_intersect = set(bins_target1.keys()).intersection(bins_target2.keys())
    if len(bin_intersect) == 0:
        return 0, None, None, None
    
    # Sort bins to ensure consistent ordering
    sorted_bins = sorted(bin_intersect)
    
    # Keep intensities of overlapping bins in consistent order
    target1_intensities = np.array([bins_target1[key] for key in sorted_bins])
    target2_intensities = np.array([bins_target2[key] for key in sorted_bins])

    # Avoid division by zero or invalid regression
    if len(target1_intensities) < 2 or np.all(target1_intensities == 0) or np.all(target2_intensities == 0):
        return 0, None, None, None

    # Fit a linear regression model
    try:
        reg = LinearRegression().fit(
            target1_intensities.reshape(-1, 1), target2_intensities
        )

        # predict y values using the model
        y_pred = reg.predict(target1_intensities.reshape(-1, 1))

        # calculate R-squared (coefficient of determination) as a measure of goodness of fit
        r2 = reg.score(target1_intensities.reshape(-1, 1), target2_intensities)

        return r2, y_pred, target1_intensities, target2_intensities
    except Exception as e:
        # Log the error and return 0
        import logging
        logger = logging.getLogger(__name__)
        logger.error(f"LinearRegression failed in binned calculation: {e}")
        return 0, None, None, None