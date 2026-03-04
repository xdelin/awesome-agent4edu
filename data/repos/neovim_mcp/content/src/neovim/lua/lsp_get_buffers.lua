local all_bufs = vim.api.nvim_list_bufs()
local ans = {}

for _, id in ipairs(all_bufs) do
    local name = vim.api.nvim_buf_get_name(id)
    local line_count = vim.api.nvim_buf_line_count(id)
    table.insert(ans, {
        id = id,
        name = name,
        line_count = line_count,
    })
end

return vim.json.encode(ans)
