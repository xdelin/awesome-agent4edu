local clients = vim.lsp.get_clients()
local result = vim.tbl_map(function(client)
    return { name = client.name, id = client.id }
end, clients)
return vim.json.encode(result)
