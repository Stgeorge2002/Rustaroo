//! Default argument configuration
//!
//!
//! Sets defaults based on mode (strict/moderate/sensitive) if they have not been
//! overwritten by the user.

use std::cmp::max;

#[derive(Debug, Clone)]
pub struct RustarooArgs {
    pub mode: String,
    pub id: Option<f64>,
    pub family_threshold: Option<f64>,
    pub len_dif_percent: Option<f64>,
    pub min_trailing_support: Option<usize>,
    pub trailing_recursive: Option<usize>,
    pub min_edge_support_sv: Option<usize>,
    pub remove_by_consensus: Option<bool>,
    pub edge_support_threshold: Option<usize>,
    pub input_files: Vec<String>,
}

impl RustarooArgs {
    pub fn new(mode: String, input_files: Vec<String>) -> Self {
        Self {
            mode,
            id: None,
            family_threshold: None,
            len_dif_percent: None,
            min_trailing_support: None,
            trailing_recursive: None,
            min_edge_support_sv: None,
            remove_by_consensus: None,
            edge_support_threshold: None,
            input_files,
        }
    }
}

/// Set defaults based on mode if they have not been overwritten by the user
pub fn set_default_args(mut args: RustarooArgs) -> RustarooArgs {
    let n_samples = args.input_files.len();

    match args.mode.as_str() {
        "strict" => {
            if args.id.is_none() {
                args.id = Some(0.98);
            }
            if args.family_threshold.is_none() {
                args.family_threshold = Some(0.7);
            }
            if args.len_dif_percent.is_none() {
                args.len_dif_percent = Some(0.98);
            }
            if args.min_trailing_support.is_none() {
                args.min_trailing_support = Some(max(2, (0.05 * n_samples as f64).ceil() as usize));
            }
            if args.trailing_recursive.is_none() {
                args.trailing_recursive = Some(99999999);
            }
            if args.min_edge_support_sv.is_none() {
                args.min_edge_support_sv = Some(max(2, (0.01 * n_samples as f64).ceil() as usize));
            }
            if args.remove_by_consensus.is_none() {
                args.remove_by_consensus = Some(true);
            }
            if args.edge_support_threshold.is_none() {
                args.edge_support_threshold =
                    Some(max(2, (0.01 * n_samples as f64).ceil() as usize));
            }
        }
        "moderate" => {
            if args.id.is_none() {
                args.id = Some(0.98);
            }
            if args.family_threshold.is_none() {
                args.family_threshold = Some(0.7);
            }
            if args.len_dif_percent.is_none() {
                args.len_dif_percent = Some(0.98);
            }
            if args.min_trailing_support.is_none() {
                args.min_trailing_support = Some(max(2, (0.01 * n_samples as f64).ceil() as usize));
            }
            if args.trailing_recursive.is_none() {
                args.trailing_recursive = Some(99999999);
            }
            if args.min_edge_support_sv.is_none() {
                args.min_edge_support_sv = Some(max(2, (0.01 * n_samples as f64).ceil() as usize));
            }
            if args.remove_by_consensus.is_none() {
                args.remove_by_consensus = Some(false);
            }
            if args.edge_support_threshold.is_none() {
                args.edge_support_threshold =
                    Some(max(2, (0.01 * n_samples as f64).ceil() as usize));
            }
        }
        "sensitive" => {
            if args.id.is_none() {
                args.id = Some(0.98);
            }
            if args.family_threshold.is_none() {
                args.family_threshold = Some(0.7);
            }
            if args.len_dif_percent.is_none() {
                args.len_dif_percent = Some(0.98);
            }
            if args.min_trailing_support.is_none() {
                args.min_trailing_support = Some(2);
            }
            if args.trailing_recursive.is_none() {
                args.trailing_recursive = Some(0);
            }
            if args.min_edge_support_sv.is_none() {
                args.min_edge_support_sv = Some(2);
            }
            if args.remove_by_consensus.is_none() {
                args.remove_by_consensus = Some(false);
            }
            if args.edge_support_threshold.is_none() {
                args.edge_support_threshold = Some(0);
            }
        }
        _ => {} // Unknown mode, don't set defaults
    }

    args
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strict_mode() {
        let args = RustarooArgs::new("strict".to_string(), vec!["file1".to_string(); 100]);
        let args = set_default_args(args);

        assert_eq!(args.id, Some(0.98));
        assert_eq!(args.family_threshold, Some(0.7));
        assert_eq!(args.len_dif_percent, Some(0.98));
        assert_eq!(args.min_trailing_support, Some(5)); // max(2, ceil(0.05 * 100))
        assert_eq!(args.trailing_recursive, Some(99999999));
        assert_eq!(args.remove_by_consensus, Some(true));
    }

    #[test]
    fn test_moderate_mode() {
        let args = RustarooArgs::new("moderate".to_string(), vec!["file1".to_string(); 100]);
        let args = set_default_args(args);

        assert_eq!(args.id, Some(0.98));
        assert_eq!(args.remove_by_consensus, Some(false));
    }

    #[test]
    fn test_sensitive_mode() {
        let args = RustarooArgs::new("sensitive".to_string(), vec!["file1".to_string(); 100]);
        let args = set_default_args(args);

        assert_eq!(args.min_trailing_support, Some(2));
        assert_eq!(args.trailing_recursive, Some(0));
        assert_eq!(args.edge_support_threshold, Some(0));
    }

    #[test]
    fn test_preserves_user_values() {
        let mut args = RustarooArgs::new("strict".to_string(), vec!["file1".to_string(); 100]);
        args.id = Some(0.95);
        args.min_trailing_support = Some(10);

        let args = set_default_args(args);

        // User values should be preserved
        assert_eq!(args.id, Some(0.95));
        assert_eq!(args.min_trailing_support, Some(10));
        // Other values should still get defaults
        assert_eq!(args.family_threshold, Some(0.7));
    }
}
