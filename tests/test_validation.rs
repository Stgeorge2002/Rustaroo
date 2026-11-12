//! Tests for validation module

use rustaroo::validation::{del_dups, escape_xml, unescape_xml};

#[test]
fn test_escape_xml() {
    assert_eq!(escape_xml("test&data"), "test&amp;data");
    assert_eq!(escape_xml("test<data>"), "test&lt;data&gt;");
    assert_eq!(escape_xml("test\"data\""), "test&quot;data&quot;");
    assert_eq!(escape_xml("test'data'"), "test&apos;data&apos;");
}

#[test]
fn test_unescape_xml() {
    assert_eq!(unescape_xml("test&amp;data"), "test&data");
    assert_eq!(unescape_xml("test&lt;data&gt;"), "test<data>");
    assert_eq!(unescape_xml("test&quot;data&quot;"), "test\"data\"");
    assert_eq!(unescape_xml("test&apos;data&apos;"), "test'data'");
}

#[test]
fn test_del_dups() {
    let input = vec!["a", "b", "a", "c", "b", "d"];
    let result = del_dups(input);
    assert_eq!(result, vec!["a", "b", "c", "d"]);
}

#[test]
fn test_del_dups_already_unique() {
    let input = vec!["a", "b", "c"];
    let result = del_dups(input);
    assert_eq!(result, vec!["a", "b", "c"]);
}

#[test]
fn test_del_dups_all_same() {
    let input = vec!["a", "a", "a"];
    let result = del_dups(input);
    assert_eq!(result, vec!["a"]);
}
