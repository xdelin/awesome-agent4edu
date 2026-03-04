import httpx
import os
import asyncio
from typing import Dict, Any, List


class WolframAlphaJSONClient:
    """Wolfram Alpha client using JSON API to avoid XML parsing issues"""
    
    def __init__(self, app_id: str):
        self.app_id = app_id
        self.base_url = "https://api.wolframalpha.com/v2/query"
    
    async def aquery(self, query: str) -> Dict[str, Any]:
        """
        Query Wolfram Alpha using JSON APId
        """
        params = {
            "input": query,
            "appid": self.app_id,
            "output": "json"
        }
        
        async with httpx.AsyncClient() as client:
            response = await client.get(self.base_url, params=params)
            response.raise_for_status()
            
            return response.json()


# Create a compatibility layer to match the original wolframalpha interface
class WolframResult:
    def __init__(self, json_result: Dict[str, Any]):
        self.json_data = json_result
        self.queryresult = json_result.get("queryresult", {})
    
    @property
    def pods(self):
        """Generator that yields Pod objects"""
        for pod_data in self.queryresult.get("pods", []):
            yield Pod(pod_data)


class Pod:
    def __init__(self, pod_data: Dict[str, Any]):
        self.data = pod_data
        self.title = pod_data.get("title", "")
    
    @property
    def subpods(self):
        """Generator that yields Subpod objects"""
        for subpod_data in self.data.get("subpods", []):
            yield Subpod(subpod_data)


class Subpod:
    def __init__(self, subpod_data: Dict[str, Any]):
        self.data = subpod_data
        self.plaintext = subpod_data.get("plaintext", "")
        
        # Handle image data if present
        img_data = subpod_data.get("img", {})
        if img_data:
            self.img = img_data
        else:
            self.img = None


class CompatibleWolframClient:
    """Client that provides compatibility with the original wolframalpha interface"""
    
    def __init__(self, app_id: str):
        self.json_client = WolframAlphaJSONClient(app_id)
    
    async def aquery(self, query: str) -> WolframResult:
        """
        Async query that returns a result compatible with the original interface
        """
        json_result = await self.json_client.aquery(query)
        return WolframResult(json_result)
    


api_key = os.getenv("WOLFRAM_API_KEY")

if api_key is None:
    raise ValueError("WOLFRAM_API_KEY environment variable not set")

client: CompatibleWolframClient = CompatibleWolframClient(api_key)

# test case for debugging your api key
if __name__ == "__main__":
    async def test():
        try:
            result = await client.aquery("1+1")
            print("Query successful!")
            
            pod_count = 0
            for pod in result.pods:
                pod_count += 1
                print(f"Pod {pod_count}: {pod.title}")
                for subpod in pod.subpods:
                    if subpod.plaintext:
                        print(f"  Result: {subpod.plaintext}")
                        
            print(f"Total pods found: {pod_count}")
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
    
    asyncio.run(test())
