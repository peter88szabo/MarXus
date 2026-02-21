//! Minimal, general-purpose text input reader for MarXus examples and small tools.
//!
//! This parser is intentionally simple:
//! - `[section]` headers start new sections (case-insensitive; stored lowercase).
//! - `key = value` pairs are stored verbatim (keys stored as given; lookups are lowercased).
//! - `#` starts a comment (ignored).
//! - Values can be continued on following lines (until next key or section).
//!
//! It is *not* intended to be a full config language; it exists so physics modules
//! can be driven by small "input decks" without embedding a huge amount of ad-hoc parsing.

use std::collections::HashMap;

pub type InputMap = HashMap<(String, String), String>;

pub fn parse_sectioned_kv(input: &str) -> Result<InputMap, String> {
    let mut map: InputMap = HashMap::new();
    let mut section = String::new();

    let mut current_key: Option<String> = None;
    let mut current_value = String::new();

    for (idx, raw_line) in input.lines().enumerate() {
        let line_no = idx + 1;

        let mut line = raw_line.trim().to_string();
        if let Some(hash) = line.find('#') {
            line = line[..hash].trim().to_string();
        }
        if line.is_empty() {
            continue;
        }

        if line.starts_with('[') && line.ends_with(']') {
            if let Some(key) = current_key.take() {
                map.insert((section.clone(), key), current_value.trim().to_string());
                current_value.clear();
            }

            section = line
                .trim_start_matches('[')
                .trim_end_matches(']')
                .trim()
                .to_lowercase();
            continue;
        }

        if let Some(eq) = line.find('=') {
            if let Some(key) = current_key.take() {
                map.insert((section.clone(), key), current_value.trim().to_string());
                current_value.clear();
            }

            let key = line[..eq].trim().to_string();
            let value = line[eq + 1..].trim().to_string();
            current_key = Some(key);
            current_value.push_str(&value);
            current_value.push('\n');
        } else if current_key.is_some() {
            current_value.push_str(&line);
            current_value.push('\n');
        } else {
            return Err(format!(
                "Invalid line (expected section header or key=value) at {}: {}",
                line_no, raw_line
            ));
        }
    }

    if let Some(key) = current_key.take() {
        map.insert((section.clone(), key), current_value.trim().to_string());
    }

    Ok(map)
}

pub fn get_string(map: &InputMap, section: &str, key: &str) -> Result<String, String> {
    let section = section.to_lowercase();
    let key = key.to_string();
    map.get(&(section.clone(), key.clone()))
        .map(|v| v.trim().to_string())
        .ok_or_else(|| format!("Missing {}.{}", section, key))
}

pub fn get_string_or(map: &InputMap, section: &str, key: &str, default: &str) -> String {
    let section = section.to_lowercase();
    let key = key.to_string();
    map.get(&(section, key))
        .map(|v| v.trim().to_string())
        .unwrap_or_else(|| default.to_string())
}

pub fn get_f64(map: &InputMap, section: &str, key: &str) -> Result<f64, String> {
    get_string(map, section, key)?
        .parse::<f64>()
        .map_err(|_| format!("Invalid f64 for {}.{}", section.to_lowercase(), key))
}

pub fn get_f64_or(map: &InputMap, section: &str, key: &str, default: f64) -> Result<f64, String> {
    let section = section.to_lowercase();
    let key = key.to_string();
    match map.get(&(section.clone(), key.clone())) {
        Some(v) => v
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("Invalid f64 for {}.{}", section, key)),
        None => Ok(default),
    }
}

pub fn get_usize(map: &InputMap, section: &str, key: &str) -> Result<usize, String> {
    get_string(map, section, key)?
        .parse::<usize>()
        .map_err(|_| format!("Invalid usize for {}.{}", section.to_lowercase(), key))
}

pub fn get_usize_or(
    map: &InputMap,
    section: &str,
    key: &str,
    default: usize,
) -> Result<usize, String> {
    let section = section.to_lowercase();
    let key = key.to_string();
    match map.get(&(section.clone(), key.clone())) {
        Some(v) => v
            .trim()
            .parse::<usize>()
            .map_err(|_| format!("Invalid usize for {}.{}", section, key)),
        None => Ok(default),
    }
}

pub fn get_vec_f64(map: &InputMap, section: &str, key: &str) -> Result<Vec<f64>, String> {
    let s = get_string(map, section, key)?;
    parse_vec_f64(&s).map_err(|e| format!("Invalid vec for {}.{}: {}", section, key, e))
}

pub fn parse_vec_f64(raw: &str) -> Result<Vec<f64>, String> {
    let s = raw.trim();
    let inner = s
        .strip_prefix('[')
        .and_then(|x| x.strip_suffix(']'))
        .ok_or_else(|| "Expected [..] list".to_string())?;
    let inner = inner.trim();
    if inner.is_empty() {
        return Ok(Vec::new());
    }
    inner
        .split(',')
        .map(|x| {
            x.trim()
                .parse::<f64>()
                .map_err(|_| format!("Invalid f64 '{}'", x.trim()))
        })
        .collect()
}
