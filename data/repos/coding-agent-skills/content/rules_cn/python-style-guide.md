---
description: 
globs: *.py,*.pyi
alwaysApply: false
---

# Python 代码规范（PEP 8 + 现代最佳实践）

## 关键规则

- **类型安全**：所有函数参数和返回值必须添加类型提示（`typing`），禁止使用 `Any`，必须启用 `mypy` 严格模式
- **异常处理**：捕获具体异常类型，禁止裸 `except:` 或 `except Exception:`，必须包含错误日志和上下文
- **性能优化**：优先使用生成器（`yield`）处理大数据，列表推导式替代简单循环，字符串拼接使用 `join`
- **不可变性**：优先使用 `final`、`frozen dataclass`，函数参数不要修改可变默认值（`def foo(x=None)`）
- **Pythonic**：善用上下文管理器（`with`）、装饰器、`@property`，避免手动管理资源（关闭文件等）
- **文档规范**：公共 API 必须包含 Google Style 或 NumPy Style 的 docstring，说明参数类型和异常

## 命名规范

- 类名：`PascalCase`（`HttpClient`, `DataProcessor`）
- 函数/变量/属性：`snake_case`（`fetch_data()`, `user_count`）
- 常量：`UPPER_SNAKE_CASE`（`MAX_RETRIES`, `DEFAULT_TIMEOUT`）
- 私有成员：单下划线前缀（`_internal_cache`），强私有使用双下划线（`__private`）
- 模块名：短横线或下划线（`data_utils.py`），包名使用下划线

## 代码结构

- 单一职责：函数不超过 20 行，类不超过 5 个公有方法，模块不超过 500 行
- 导入排序：标准库 → 第三方库 → 本地模块，每组空行分隔，使用 `isort` 自动化
- 避免循环导入：使用 `TYPE_CHECKING` 进行类型检查时的导入，或重构模块结构
- 严格相等：使用 `is` 判断 `None`（`if x is None`），禁止 `== None`

## 类型系统

- 强制类型提示：`def process(data: dict[str, int]) -&gt; list[str]:`
- 复杂类型使用别名：`UserDict = dict[str, UserInfo]`
- 可选参数明确标注：`def find(x: int | None = None) -&gt; User | None:`
- 泛型约束：`T = TypeVar('T', bound=BaseModel)`

## 安全与鲁棒性

- 输入验证：使用 `pydantic` 或 `marshmallow` 验证外部数据，禁止直接解包 `**kwargs`
- 路径安全：使用 `pathlib.Path` 替代字符串拼接路径，禁止 `os.path.join`
- SQL 安全：必须使用参数化查询（SQLAlchemy `text()` + bind params），禁止 f-string 拼接 SQL
- 环境配置：使用 `pydantic-settings` 或 `python-dotenv`，禁止硬编码密钥

&lt;example&gt;
# 正确：类型提示、异常处理、生成器、Pythonic
from pathlib import Path
from typing import Iterator, Final
import logging
import json

logger = logging.getLogger(__name__)

CHUNK_SIZE: Final[int] = 8192

def process_large_file(file_path: Path) -&gt; Iterator[dict]:
    """
    逐行处理大 JSON 文件，返回生成器避免内存爆炸。
    
    Args:
        file_path: JSON 文件路径，每行一个 JSON 对象
        
    Yields:
        解析后的字典对象
        
    Raises:
        FileNotFoundError: 文件不存在
        json.JSONDecodeError: JSON 格式错误
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    try:
        with file_path.open('r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError as e:
                    logger.error(f"JSON decode error at line {line_num}: {e}")
                    raise
    except OSError as e:
        logger.error(f"File operation failed: {e}")
        raise RuntimeError(f"Cannot process file {file_path}") from e

# 使用：内存友好
for record in process_large_file(Path("data.jsonl")):
    process_record(record)
&lt;/example&gt;

&lt;example type="invalid"&gt;
# 错误：无类型提示、裸异常、内存爆炸、字符串路径
import os
import json

def process_file(path):
    # 错误：未验证路径，使用字符串而非 Path
    f = open(os.path.join("/data", path), "r")  # 资源泄漏风险
    
    lines = f.readlines()  # 内存爆炸：一次性读入所有行
    
    results = []
    for i in range(len(lines)):  # 非 Pythonic，应该用 enumerate
        try:
            data = json.loads(lines[i])
            results.append(data)
        except:  # 裸异常，捕获 KeyboardInterrupt 等
            pass  # 静默失败
    
    f.close()  # 如果上面异常，这里不执行，句柄泄漏
    return results
&lt;/example&gt;