//! Tests for argument handling

use rustaroo::args::{set_default_args, RustarooArgs};

#[test]
fn test_default_args_strict() {
    let args = RustarooArgs::new("strict".to_string(), vec![]);
    let args = set_default_args(args);

    assert!(args.id.is_some());
    assert!(args.family_threshold.is_some());
    assert!(args.len_dif_percent.is_some());

    // Strict mode should have higher thresholds
    assert!(args.id.unwrap() >= 0.95);
}

#[test]
fn test_default_args_moderate() {
    let args = RustarooArgs::new("moderate".to_string(), vec![]);
    let args = set_default_args(args);

    assert!(args.id.is_some());
    assert!(args.family_threshold.is_some());
}

#[test]
fn test_default_args_sensitive() {
    let args = RustarooArgs::new("sensitive".to_string(), vec![]);
    let args = set_default_args(args);

    assert!(args.id.is_some());
    assert!(args.family_threshold.is_some());

    // Sensitive mode should have lower thresholds
    assert!(args.id.unwrap() <= 0.98);
}

#[test]
fn test_args_preserve_user_values() {
    let mut args = RustarooArgs::new("strict".to_string(), vec![]);
    args.id = Some(0.99);
    let args = set_default_args(args);

    // User-specified value should be preserved
    assert_eq!(args.id.unwrap(), 0.99);
}
