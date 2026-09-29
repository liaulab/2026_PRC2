import json, re

def _scan_value(s, i):
    """Scan one JSON value (array or object) starting at s[i], return end index (exclusive)."""
    assert s[i] in '[{'
    open_c = s[i]
    close_c = ']' if open_c == '[' else '}'
    depth = 0
    in_str = False
    esc = False
    j = i
    while j < len(s):
        c = s[j]
        if in_str:
            if esc:
                esc = False
            elif c == '\\':
                esc = True
            elif c == '"':
                in_str = False
        else:
            if c == '"':
                in_str = True
            elif c == open_c:
                depth += 1
            elif c == close_c:
                depth -= 1
                if depth == 0:
                    return j + 1
            elif c in '[{':
                # nested other-type bracket; just track generically via stack instead
                pass
        j += 1
    raise ValueError("unterminated value")

def _scan_any(s, i):
    """More general: track a stack of brackets so [ and { nesting both work."""
    stack = []
    in_str = False
    esc = False
    j = i
    start = i
    while j < len(s):
        c = s[j]
        if in_str:
            if esc:
                esc = False
            elif c == '\\':
                esc = True
            elif c == '"':
                in_str = False
        else:
            if c == '"':
                in_str = True
            elif c in '[{':
                stack.append(c)
            elif c in ']}':
                stack.pop()
                if not stack:
                    return j + 1
        j += 1
    raise ValueError("unterminated value")

def find_newplot(html):
    m = re.search(r'Plotly\.newPlot\(\s*"([^"]+)"\s*,\s*', html)
    if not m:
        raise ValueError("no Plotly.newPlot found")
    div_id = m.group(1)
    i = m.end()
    assert html[i] == '['
    end_data = _scan_any(html, i)
    data_s = html[i:end_data]
    j = end_data
    while html[j] in ' \t\n,':
        j += 1
    assert html[j] == '{'
    end_layout = _scan_any(html, j)
    layout_s = html[j:end_layout]
    k = end_layout
    while html[k] in ' \t\n,':
        k += 1
    assert html[k] == '{'
    end_config = _scan_any(html, k)
    config_s = html[k:end_config]
    # find_newplot only located the opening "Plotly.newPlot(" -- the call's
    # own closing ")" still follows end_config in the source and must be
    # consumed too, or a replacement built from these pieces (which supplies
    # its own closing paren) leaves the original ")" behind as a stray
    # token, breaking the script silently (no console-visible syntax
    # checking happens until the browser actually parses it).
    n = end_config
    while html[n] in ' \t\n':
        n += 1
    assert html[n] == ')', f"expected ')' at {n}, found {html[n:n+20]!r}"
    span = (m.start(), n + 1)
    return div_id, data_s, layout_s, config_s, span

def load(path):
    html = open(path, encoding='utf-8').read()
    div_id, data_s, layout_s, config_s, span = find_newplot(html)
    data = json.loads(data_s)
    layout = json.loads(layout_s)
    config = json.loads(config_s)
    return html, div_id, data, layout, config, span
