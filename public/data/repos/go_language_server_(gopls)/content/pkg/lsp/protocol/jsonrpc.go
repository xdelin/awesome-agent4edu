package protocol

import (
	"encoding/json"
	"fmt"
)

type JSONRPCMessage struct {
	JSONRPC string          `json:"jsonrpc"`
	ID      any             `json:"id,omitempty"`
	Method  string          `json:"method,omitempty"`
	Params  json.RawMessage `json:"params,omitempty"`
	Result  json.RawMessage `json:"result,omitempty"`
	Error   *JSONRPCError   `json:"error,omitempty"`
}

type JSONRPCError struct {
	Code    int             `json:"code"`
	Message string          `json:"message"`
	Data    json.RawMessage `json:"data,omitempty"`
}

func (e *JSONRPCError) Error() string {
	return fmt.Sprintf("JSON-RPC error %d: %s", e.Code, e.Message)
}

func (m *JSONRPCMessage) UnmarshalJSON(data []byte) error {
	type Alias JSONRPCMessage
	aux := struct {
		ID json.RawMessage `json:"id,omitempty"`
		*Alias
	}{
		Alias: (*Alias)(m),
	}

	if err := json.Unmarshal(data, &aux); err != nil {
		return err
	}

	if len(aux.ID) > 0 {
		var id any
		var num json.Number
		if err := json.Unmarshal(aux.ID, &num); err == nil {
			id = num
		} else {
			var str string
			if json.Unmarshal(aux.ID, &str) == nil {
				id = str
			} else {
				id = aux.ID
			}
		}
		m.ID = id
	}

	return nil
}

func NewRequest(id any, method string, params any) (*JSONRPCMessage, error) {
	var paramsRaw json.RawMessage
	if params != nil {
		var err error
		paramsRaw, err = json.Marshal(params)
		if err != nil {
			return nil, fmt.Errorf("failed to marshal params: %w", err)
		}
	}

	return &JSONRPCMessage{
		JSONRPC: "2.0",
		ID:      id,
		Method:  method,
		Params:  paramsRaw,
	}, nil
}

func NewNotification(method string, params any) (*JSONRPCMessage, error) {
	var paramsRaw json.RawMessage
	if params != nil {
		var err error
		paramsRaw, err = json.Marshal(params)
		if err != nil {
			return nil, fmt.Errorf("failed to marshal params: %w", err)
		}
	}

	return &JSONRPCMessage{
		JSONRPC: "2.0",
		Method:  method,
		Params:  paramsRaw,
	}, nil
}

func (msg *JSONRPCMessage) ParseResult(target any) error {
	if msg.Error != nil {
		return msg.Error
	}

	if msg.Result == nil {
		return fmt.Errorf("no result in response")
	}

	return json.Unmarshal(msg.Result, target)
}
