---
description: 
globs: *.cpp,*.hpp,*.h,*.cc,*.cxx
alwaysApply: false
---

# C++ 代码规范（Modern C++17/20/23）

## 关键规则

- **内存安全**：禁止裸指针管理资源，必须使用智能指针（`unique_ptr`/`shared_ptr`）或容器
- **RAII**：资源获取即初始化，确保异常安全，构造函数不抛异常时使用 `noexcept`
- **现代特性**：优先使用 `auto`、`range-based for`、`constexpr`、`nullptr`，禁用 `NULL` 和 `0`
- **异常安全**：所有函数必须指定异常规格（`noexcept` 或异常说明），析构函数必须 `noexcept`
- **类型安全**：禁用 C 风格类型转换，使用 `static_cast`/`dynamic_cast`，避免 `reinterpret_cast`
- **模板约束**：模板参数必须使用 `concept`（C++20）或 `static_assert` 约束，提供清晰编译错误
- **性能优化**：优先使用移动语义（`&&`）和完美转发，避免不必要的拷贝，热路径禁用虚函数
- **编译安全**：头文件必须包含 `#pragma once` 或 include guard，前向声明优于包含头文件

## 命名规范

- 类/结构体：`PascalCase`（`HttpRequest`, `JsonParser`）
- 变量/函数：`snake_case`（`process_data()`, `max_buffer_size`）
- 私有成员：尾部下划线（`private_member_`）或 `m_` 前缀（`m_private_member`），团队统一即可
- 宏/常量：`UPPER_SNAKE_CASE`（`BUFFER_SIZE`, `DEBUG_MODE`）
- 模板参数：`CamelCase`（`InputIt`, `Alloc`）

## 代码组织

- 头文件只暴露接口，实现细节放在 `.cpp`，模板类除外
- 单一职责：类不超过 5 个公有方法，函数不超过 30 行（简单委托除外）
- 包含顺序：相关头文件 → C 库 → C++ 库 → 第三方库 → 项目头文件，每类空行分隔
- 禁用 `using namespace std;`，可使用 `using std::string` 具体引入

## 内存与资源

- 裸 `new`/`delete` 零容忍，使用 `std::make_unique` 和 `std::make_shared`
- 数组使用 `std::vector` 或 `std::array`，禁止 C 风格数组（`int arr[10]`）
- 字符串使用 `std::string`，禁止 `char*` 操作（`strcpy`, `sprintf` 等）
- 锁管理必须使用 `std::lock_guard` 或 `std::unique_lock`，禁止裸 `mutex.lock()`

&lt;example&gt;
// 正确：RAII、智能指针、异常安全、现代 C++
class FileProcessor {
public:
    explicit FileProcessor(std::string_view path) 
        : file_(std::fopen(path.data(), "r")) {
        if (!file_) {
            throw std::runtime_error("Failed to open file: " + std::string(path));
        }
    }
    
    // 禁用拷贝，支持移动
    FileProcessor(const FileProcessor&) = delete;
    FileProcessor& operator=(const FileProcessor&) = delete;
    
    FileProcessor(FileProcessor&& other) noexcept 
        : file_(std::exchange(other.file_, nullptr)) {}
    
    ~FileProcessor() noexcept {
        if (file_) std::fclose(file_);
    }
    
    std::string read_all() {
        std::string content;
        // 使用现代 C++ 特性
        constexpr size_t buffer_size = 4096;
        std::vector&lt;char&gt; buffer(buffer_size);
        
        while (std::fgets(buffer.data(), buffer.size(), file_)) {
            content.append(buffer.data());
        }
        return content;
    }

private:
    std::FILE* file_;
};

// 使用
void process() {
    auto processor = std::make_unique&lt;FileProcessor&gt;("data.txt");
    auto content = processor-&gt;read_all(); // 返回值优化（RVO）
}
&lt;/example&gt;

&lt;example type="invalid"&gt;
// 错误：裸指针、内存泄漏、异常不安全、C 风格
class BadProcessor {
public:
    BadProcessor(const char* path) {
        file = fopen(path, "r"); // 可能失败未检查
        buffer = new char[1024]; // 裸 new
    }
    
    ~BadProcessor() {
        fclose(file);
        delete[] buffer; // 如果构造函数异常，这里不执行，内存泄漏
    }
    
    char* read() {
        fgets(buffer, 1024, file);
        return buffer; // 悬挂指针风险
    }
    
private:
    FILE* file;
    char* buffer;
};
&lt;/example&gt;