//! Validation utilities and GML string handling
//!

use anyhow::{bail, Result};
use regex::Regex;
use std::collections::HashSet;
use std::path::Path;

/// XML character reference escaping for GML output
pub fn escape(text: &str) -> String {
    let mut result = String::new();
    for ch in text.chars() {
        if ch.is_ascii() && ch >= ' ' && ch <= '~' && ch != '&' && ch != '"' {
            result.push(ch);
        } else {
            result.push_str(&format!("&#{};", ch as u32));
        }
    }
    result
}

/// Replace XML character references with the referenced characters
pub fn unescape(text: &str) -> String {
    let re = Regex::new(r"&(?:[0-9A-Za-z]+|#(?:[0-9]+|x[0-9A-Fa-f]+));").unwrap();

    re.replace_all(text, |caps: &regex::Captures| {
        let entity = &caps[0];
        if entity.starts_with("&#") {
            let num_part = &entity[2..entity.len() - 1];
            if num_part.starts_with('x') {
                // Hex
                if let Ok(code) = u32::from_str_radix(&num_part[1..], 16) {
                    if let Some(ch) = char::from_u32(code) {
                        return ch.to_string();
                    }
                }
            } else {
                // Decimal
                if let Ok(code) = num_part.parse::<u32>() {
                    if let Some(ch) = char::from_u32(code) {
                        return ch.to_string();
                    }
                }
            }
        } else {
            // Named entity - would need full HTML entity table
            // For now, just handle common ones
            match &entity[1..entity.len() - 1] {
                "amp" => return "&".to_string(),
                "lt" => return "<".to_string(),
                "gt" => return ">".to_string(),
                "quot" => return "\"".to_string(),
                "apos" => return "'".to_string(),
                _ => {}
            }
        }
        entity.to_string() // leave unchanged if can't decode
    })
    .to_string()
}

/// Escape text using standard XML/named character references.
pub fn escape_xml(text: &str) -> String {
    let mut result = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '&' => result.push_str("&amp;"),
            '<' => result.push_str("&lt;"),
            '>' => result.push_str("&gt;"),
            '"' => result.push_str("&quot;"),
            '\'' => result.push_str("&apos;"),
            _ => result.push(ch),
        }
    }
    result
}

/// Replace standard XML/named and numeric character references.
pub fn unescape_xml(text: &str) -> String {
    unescape(text)
}

/// Validate that a file exists
pub fn is_valid_file(path: &str) -> Result<String> {
    if !Path::new(path).exists() {
        bail!("The file {} does not exist!", path);
    }
    Ok(path.to_string())
}

/// Validate that a folder exists
pub fn is_valid_folder(path: &str) -> Result<String> {
    if !Path::new(path).is_dir() {
        bail!("The folder {} does not exist!", path);
    }
    Ok(path.to_string())
}

/// Convert a value to list if it isn't already
pub fn conv_list<T>(maybe_list: Vec<T>) -> Vec<T> {
    maybe_list
}

/// Remove duplicates from a sequence
pub fn del_dups<T: Eq + std::hash::Hash + Clone>(seq: Vec<T>) -> Vec<T> {
    let mut seen = HashSet::new();
    let mut result = Vec::new();
    for item in seq {
        if seen.insert(item.clone()) {
            result.push(item);
        }
    }
    result
}

/// Custom stringizer for GML representation
///
/// Converts values to Python literal strings in GML format,
/// matching the behavior of the original Python implementation.
pub fn custom_stringizer(value: &serde_json::Value) -> Result<String> {
    fn stringize(value: &serde_json::Value, buf: &mut String) -> Result<()> {
        match value {
            serde_json::Value::Null => buf.push_str("None"),
            serde_json::Value::Bool(true) => buf.push('1'),
            serde_json::Value::Bool(false) => buf.push('0'),
            serde_json::Value::Number(n) => {
                if n.is_i64() {
                    buf.push_str(&n.as_i64().unwrap().to_string());
                } else if n.is_f64() {
                    buf.push_str(&n.as_f64().unwrap().to_string());
                } else {
                    buf.push_str(&n.to_string());
                }
            }
            serde_json::Value::String(s) => {
                buf.push_str(&format!("{:?}", s));
            }
            serde_json::Value::Array(arr) => {
                buf.push('[');
                for (i, item) in arr.iter().enumerate() {
                    if i > 0 {
                        buf.push(',');
                    }
                    stringize(item, buf)?;
                }
                buf.push(']');
            }
            serde_json::Value::Object(obj) => {
                buf.push('{');
                for (i, (key, val)) in obj.iter().enumerate() {
                    if i > 0 {
                        buf.push(',');
                    }
                    buf.push_str(&format!("{:?}", key));
                    buf.push(':');
                    stringize(val, buf)?;
                }
                buf.push('}');
            }
        }
        Ok(())
    }

    let mut buf = String::new();
    stringize(value, &mut buf)?;
    Ok(buf)
}

/// Check if a gene sequence is valid
///
/// Validates that:
/// - DNA sequence length is divisible by 3
/// - Protein starts with 'M' (after stripping 'X')
/// - No premature stop codons in protein sequence
pub fn is_valid_gene(dna: &str, protein: &str) -> bool {
    // Check if sequence is divisible by 3
    if dna.len() % 3 != 0 {
        return false;
    }

    // Check for start codon
    let protein_trimmed = protein.trim_matches('X');
    if protein_trimmed.is_empty() || !protein_trimmed.starts_with('M') {
        return false;
    }

    // Check for premature stop
    if let Some(last_char) = protein_trimmed.chars().last() {
        let without_last = &protein_trimmed[..protein_trimmed.len() - last_char.len_utf8()];
        if without_last.contains('*') {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_escape() {
        assert_eq!(escape("hello"), "hello");
        assert_eq!(escape("hello world"), "hello world");
        assert_eq!(escape("test&test"), "test&#38;test");
        assert_eq!(escape("test\"test"), "test&#34;test");
    }

    #[test]
    fn test_unescape() {
        assert_eq!(unescape("hello"), "hello");
        assert_eq!(unescape("&#38;"), "&");
        assert_eq!(unescape("&#34;"), "\"");
        assert_eq!(unescape("&amp;"), "&");
    }

    #[test]
    fn test_is_valid_gene() {
        // Valid gene
        assert!(is_valid_gene("ATGTAG", "M*"));
        assert!(is_valid_gene("ATGGCCTAG", "MA*"));

        // Invalid: not divisible by 3
        assert!(!is_valid_gene("ATGTA", "M"));

        // Invalid: doesn't start with M
        assert!(!is_valid_gene("GCCTAG", "A*"));

        // Invalid: premature stop
        assert!(!is_valid_gene("ATGTAGGCC", "M*A"));
    }

    #[test]
    fn test_del_dups() {
        assert_eq!(del_dups(vec![1, 2, 2, 3, 1, 4]), vec![1, 2, 3, 4]);
        assert_eq!(del_dups(vec!["a", "b", "a"]), vec!["a", "b"]);
    }
}
