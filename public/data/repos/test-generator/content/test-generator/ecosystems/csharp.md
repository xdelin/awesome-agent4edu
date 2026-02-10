# C# / .NET Testing Ecosystem

Comprehensive reference for testing C# and .NET projects.

## Detection

**Manifest files:** `*.csproj`, `*.sln`, `Directory.Build.props`
**Test frameworks:** xUnit, NUnit, MSTest

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `xunit` in csproj | xUnit |
| `NUnit` in csproj | NUnit |
| `MSTest.TestFramework` in csproj | MSTest |
| `Moq` in csproj | Moq (mocking) |
| `NSubstitute` in csproj | NSubstitute (mocking) |
| `FluentAssertions` in csproj | FluentAssertions |

---

## File Structure

### Naming Conventions
```
*Tests.cs           # Standard suffix
*Test.cs            # Alternative
*Specs.cs           # BDD style
```

### Project Structure

**Separate test project (standard):**
```
Solution/
├── src/
│   └── MyApp/
│       ├── MyApp.csproj
│       ├── Services/
│       │   └── UserService.cs
│       └── Models/
│           └── User.cs
└── tests/
    └── MyApp.Tests/
        ├── MyApp.Tests.csproj
        ├── Services/
        │   └── UserServiceTests.cs
        └── Models/
            └── UserTests.cs
```

**With integration tests:**
```
tests/
├── MyApp.UnitTests/
│   └── MyApp.UnitTests.csproj
├── MyApp.IntegrationTests/
│   └── MyApp.IntegrationTests.csproj
└── MyApp.FunctionalTests/
    └── MyApp.FunctionalTests.csproj
```

---

## xUnit Patterns

### Basic Structure
```csharp
using Xunit;

namespace MyApp.Tests.Services;

public class UserServiceTests
{
    private readonly UserService _sut;  // System Under Test

    public UserServiceTests()
    {
        // Runs before each test (constructor = setup)
        _sut = new UserService();
    }

    [Fact]
    public void GetById_ValidId_ReturnsUser()
    {
        // Arrange
        var userId = 1;

        // Act
        var result = _sut.GetById(userId);

        // Assert
        Assert.NotNull(result);
        Assert.Equal(userId, result.Id);
    }

    [Fact]
    public void GetById_InvalidId_ThrowsException()
    {
        // Arrange
        var invalidId = -1;

        // Act & Assert
        Assert.Throws<ArgumentException>(() => _sut.GetById(invalidId));
    }
}
```

### Assertions
```csharp
// Equality
Assert.Equal(expected, actual);
Assert.NotEqual(unexpected, actual);
Assert.Same(expected, actual);          // Reference equality
Assert.NotSame(unexpected, actual);

// Boolean
Assert.True(condition);
Assert.False(condition);
Assert.Null(value);
Assert.NotNull(value);

// Strings
Assert.Contains("substring", actualString);
Assert.DoesNotContain("substring", actualString);
Assert.StartsWith("prefix", actualString);
Assert.EndsWith("suffix", actualString);
Assert.Matches("regex", actualString);
Assert.Empty(actualString);

// Collections
Assert.Empty(collection);
Assert.NotEmpty(collection);
Assert.Contains(item, collection);
Assert.DoesNotContain(item, collection);
Assert.Single(collection);
Assert.All(collection, item => Assert.True(item > 0));

// Types
Assert.IsType<ExpectedType>(actual);
Assert.IsAssignableFrom<BaseType>(actual);

// Exceptions
var ex = Assert.Throws<ArgumentException>(() => Method());
Assert.Equal("Expected message", ex.Message);

await Assert.ThrowsAsync<InvalidOperationException>(
    async () => await AsyncMethod());

// Ranges
Assert.InRange(actual, low, high);
```

### Theory (Parameterized Tests)
```csharp
public class CalculatorTests
{
    [Theory]
    [InlineData(1, 2, 3)]
    [InlineData(0, 0, 0)]
    [InlineData(-1, 1, 0)]
    [InlineData(100, 200, 300)]
    public void Add_TwoNumbers_ReturnsSum(int a, int b, int expected)
    {
        var calculator = new Calculator();
        var result = calculator.Add(a, b);
        Assert.Equal(expected, result);
    }

    [Theory]
    [MemberData(nameof(GetValidationTestData))]
    public void Validate_TestCases(string input, bool expected)
    {
        var result = Validator.IsValid(input);
        Assert.Equal(expected, result);
    }

    public static IEnumerable<object[]> GetValidationTestData()
    {
        yield return new object[] { "valid@email.com", true };
        yield return new object[] { "invalid-email", false };
        yield return new object[] { "", false };
        yield return new object[] { null, false };
    }

    [Theory]
    [ClassData(typeof(CalculatorTestData))]
    public void Add_FromClassData_ReturnsSum(int a, int b, int expected)
    {
        // Uses data from separate class
    }
}

public class CalculatorTestData : IEnumerable<object[]>
{
    public IEnumerator<object[]> GetEnumerator()
    {
        yield return new object[] { 1, 2, 3 };
        yield return new object[] { -1, -1, -2 };
    }

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}
```

### Test Fixtures (Shared Context)
```csharp
public class DatabaseFixture : IDisposable
{
    public DbConnection Connection { get; }

    public DatabaseFixture()
    {
        Connection = new SqlConnection("...");
        Connection.Open();
        // Seed data
    }

    public void Dispose()
    {
        Connection.Close();
    }
}

[CollectionDefinition("Database collection")]
public class DatabaseCollection : ICollectionFixture<DatabaseFixture>
{
}

[Collection("Database collection")]
public class DatabaseTests
{
    private readonly DatabaseFixture _fixture;

    public DatabaseTests(DatabaseFixture fixture)
    {
        _fixture = fixture;
    }

    [Fact]
    public void Query_ReturnsData()
    {
        // Use _fixture.Connection
    }
}
```

---

## NUnit Patterns

### Basic Structure
```csharp
using NUnit.Framework;

namespace MyApp.Tests.Services;

[TestFixture]
public class UserServiceTests
{
    private UserService _sut;

    [SetUp]
    public void Setup()
    {
        _sut = new UserService();
    }

    [TearDown]
    public void TearDown()
    {
        // Cleanup
    }

    [OneTimeSetUp]
    public void OneTimeSetUp()
    {
        // Run once before all tests
    }

    [OneTimeTearDown]
    public void OneTimeTearDown()
    {
        // Run once after all tests
    }

    [Test]
    public void GetById_ValidId_ReturnsUser()
    {
        // Arrange
        var userId = 1;

        // Act
        var result = _sut.GetById(userId);

        // Assert
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Id, Is.EqualTo(userId));
    }
}
```

### NUnit Assertions (Constraint Model)
```csharp
// Equality
Assert.That(actual, Is.EqualTo(expected));
Assert.That(actual, Is.Not.EqualTo(unexpected));
Assert.That(actual, Is.SameAs(expected));

// Boolean
Assert.That(condition, Is.True);
Assert.That(condition, Is.False);
Assert.That(value, Is.Null);
Assert.That(value, Is.Not.Null);

// Numbers
Assert.That(value, Is.GreaterThan(5));
Assert.That(value, Is.LessThanOrEqualTo(10));
Assert.That(value, Is.InRange(1, 100));
Assert.That(value, Is.Positive);
Assert.That(value, Is.Negative);

// Strings
Assert.That(str, Does.Contain("substring"));
Assert.That(str, Does.StartWith("prefix"));
Assert.That(str, Does.EndWith("suffix"));
Assert.That(str, Does.Match("regex"));
Assert.That(str, Is.Empty);

// Collections
Assert.That(collection, Is.Empty);
Assert.That(collection, Is.Not.Empty);
Assert.That(collection, Has.Count.EqualTo(3));
Assert.That(collection, Does.Contain(item));
Assert.That(collection, Is.All.GreaterThan(0));
Assert.That(collection, Is.Unique);
Assert.That(collection, Is.Ordered);

// Types
Assert.That(actual, Is.InstanceOf<ExpectedType>());
Assert.That(actual, Is.AssignableTo<BaseType>());

// Exceptions
Assert.That(() => Method(), Throws.TypeOf<ArgumentException>());
Assert.That(() => Method(), Throws.ArgumentException);
Assert.That(() => Method(), 
    Throws.TypeOf<ArgumentException>()
        .With.Message.Contains("invalid"));

// Async
Assert.That(async () => await AsyncMethod(), 
    Throws.TypeOf<InvalidOperationException>());
```

### NUnit Parameterized Tests
```csharp
[TestFixture]
public class CalculatorTests
{
    [TestCase(1, 2, 3)]
    [TestCase(0, 0, 0)]
    [TestCase(-1, 1, 0)]
    public void Add_TwoNumbers_ReturnsSum(int a, int b, int expected)
    {
        var result = Calculator.Add(a, b);
        Assert.That(result, Is.EqualTo(expected));
    }

    [TestCaseSource(nameof(DivisionTestCases))]
    public void Divide_TestCases(int a, int b, double expected)
    {
        var result = Calculator.Divide(a, b);
        Assert.That(result, Is.EqualTo(expected).Within(0.001));
    }

    private static IEnumerable<TestCaseData> DivisionTestCases()
    {
        yield return new TestCaseData(10, 2, 5.0).SetName("Divide 10 by 2");
        yield return new TestCaseData(7, 2, 3.5).SetName("Divide 7 by 2");
    }

    [Test]
    [TestCase("", ExpectedResult = false)]
    [TestCase("valid@email.com", ExpectedResult = true)]
    [TestCase("invalid", ExpectedResult = false)]
    public bool IsValidEmail_TestCases(string email)
    {
        return Validator.IsValidEmail(email);
    }
}
```

---

## MSTest Patterns

### Basic Structure
```csharp
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MyApp.Tests.Services;

[TestClass]
public class UserServiceTests
{
    private UserService _sut;

    [TestInitialize]
    public void Setup()
    {
        _sut = new UserService();
    }

    [TestCleanup]
    public void Cleanup()
    {
        // Cleanup
    }

    [ClassInitialize]
    public static void ClassSetup(TestContext context)
    {
        // Run once before all tests
    }

    [ClassCleanup]
    public static void ClassCleanup()
    {
        // Run once after all tests
    }

    [TestMethod]
    public void GetById_ValidId_ReturnsUser()
    {
        var result = _sut.GetById(1);
        
        Assert.IsNotNull(result);
        Assert.AreEqual(1, result.Id);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void GetById_InvalidId_ThrowsException()
    {
        _sut.GetById(-1);
    }
}
```

### MSTest Data-Driven Tests
```csharp
[TestClass]
public class CalculatorTests
{
    [DataTestMethod]
    [DataRow(1, 2, 3)]
    [DataRow(0, 0, 0)]
    [DataRow(-1, 1, 0)]
    public void Add_TwoNumbers_ReturnsSum(int a, int b, int expected)
    {
        var result = Calculator.Add(a, b);
        Assert.AreEqual(expected, result);
    }

    [DataTestMethod]
    [DynamicData(nameof(GetTestData), DynamicDataSourceType.Method)]
    public void Validate_TestCases(string input, bool expected)
    {
        var result = Validator.IsValid(input);
        Assert.AreEqual(expected, result);
    }

    private static IEnumerable<object[]> GetTestData()
    {
        yield return new object[] { "valid", true };
        yield return new object[] { "", false };
    }
}
```

---

## Moq (Mocking)

### Basic Mocking
```csharp
using Moq;

public class UserServiceTests
{
    private readonly Mock<IUserRepository> _mockRepository;
    private readonly Mock<IEmailService> _mockEmailService;
    private readonly UserService _sut;

    public UserServiceTests()
    {
        _mockRepository = new Mock<IUserRepository>();
        _mockEmailService = new Mock<IEmailService>();
        _sut = new UserService(_mockRepository.Object, _mockEmailService.Object);
    }

    [Fact]
    public void GetById_ValidId_ReturnsUser()
    {
        // Arrange
        var expectedUser = new User { Id = 1, Name = "John" };
        _mockRepository
            .Setup(r => r.GetById(1))
            .Returns(expectedUser);

        // Act
        var result = _sut.GetById(1);

        // Assert
        Assert.Equal("John", result.Name);
        _mockRepository.Verify(r => r.GetById(1), Times.Once);
    }

    [Fact]
    public async Task CreateUser_ValidUser_SendsEmail()
    {
        // Arrange
        var user = new User { Email = "john@example.com" };
        _mockRepository
            .Setup(r => r.AddAsync(It.IsAny<User>()))
            .ReturnsAsync(user);

        // Act
        await _sut.CreateAsync(user);

        // Assert
        _mockEmailService.Verify(
            e => e.SendWelcomeEmailAsync(user.Email),
            Times.Once);
    }
}
```

### Setup Patterns
```csharp
// Return value
mock.Setup(x => x.Method()).Returns(value);
mock.Setup(x => x.Method()).ReturnsAsync(value);  // Async

// Sequence
mock.SetupSequence(x => x.Method())
    .Returns(first)
    .Returns(second)
    .Throws<Exception>();

// Throw exception
mock.Setup(x => x.Method()).Throws<InvalidOperationException>();
mock.Setup(x => x.Method()).ThrowsAsync(new Exception());

// Callback
mock.Setup(x => x.Method(It.IsAny<string>()))
    .Callback<string>(s => Console.WriteLine(s))
    .Returns(true);

// Properties
mock.Setup(x => x.Property).Returns(value);
mock.SetupProperty(x => x.Property, defaultValue);
mock.SetupAllProperties();  // All properties trackable
```

### Argument Matching
```csharp
// Any value
mock.Setup(x => x.Method(It.IsAny<int>())).Returns(true);
mock.Setup(x => x.Method(It.IsAny<string>())).Returns(true);

// Specific conditions
mock.Setup(x => x.Method(It.Is<int>(i => i > 0))).Returns(true);
mock.Setup(x => x.Method(It.IsIn(1, 2, 3))).Returns(true);
mock.Setup(x => x.Method(It.IsInRange(1, 10, Range.Inclusive))).Returns(true);
mock.Setup(x => x.Method(It.IsRegex("pattern"))).Returns(true);

// Null
mock.Setup(x => x.Method(null)).Returns(false);
mock.Setup(x => x.Method(It.IsNotNull<string>())).Returns(true);
```

### Verification
```csharp
// Basic
mock.Verify(x => x.Method());
mock.Verify(x => x.Method(), Times.Once);
mock.Verify(x => x.Method(), Times.Exactly(2));
mock.Verify(x => x.Method(), Times.Never);
mock.Verify(x => x.Method(), Times.AtLeastOnce);
mock.Verify(x => x.Method(), Times.AtMost(3));

// With arguments
mock.Verify(x => x.Method("expected"), Times.Once);
mock.Verify(x => x.Method(It.Is<string>(s => s.Contains("test"))));

// Verify no other calls
mock.VerifyNoOtherCalls();
```

---

## FluentAssertions

```csharp
using FluentAssertions;

// Strings
name.Should().Be("John");
name.Should().StartWith("Jo").And.EndWith("hn");
name.Should().Contain("oh");
name.Should().BeNullOrEmpty();
name.Should().MatchRegex("pattern");

// Numbers
age.Should().Be(25);
age.Should().BeGreaterThan(18);
age.Should().BeInRange(18, 65);
price.Should().BeApproximately(10.0, 0.01);

// Collections
list.Should().HaveCount(3);
list.Should().Contain(item);
list.Should().ContainInOrder(1, 2, 3);
list.Should().OnlyContain(x => x > 0);
list.Should().BeInAscendingOrder();

// Objects
user.Should().BeEquivalentTo(expectedUser);
user.Should().BeEquivalentTo(expectedUser, options => 
    options.Excluding(u => u.Id));

// Exceptions
action.Should().Throw<ArgumentException>()
    .WithMessage("*invalid*");

await asyncAction.Should().ThrowAsync<InvalidOperationException>();

action.Should().NotThrow();

// Execution time
action.ExecutionTime().Should().BeLessThan(1.Seconds());
```

---

## ASP.NET Core Integration Tests

```csharp
using Microsoft.AspNetCore.Mvc.Testing;

public class UserApiTests : IClassFixture<WebApplicationFactory<Program>>
{
    private readonly HttpClient _client;

    public UserApiTests(WebApplicationFactory<Program> factory)
    {
        _client = factory.CreateClient();
    }

    [Fact]
    public async Task GetUsers_ReturnsSuccessStatusCode()
    {
        var response = await _client.GetAsync("/api/users");

        response.EnsureSuccessStatusCode();
        response.Content.Headers.ContentType?.ToString()
            .Should().Be("application/json; charset=utf-8");
    }

    [Fact]
    public async Task CreateUser_ValidUser_ReturnsCreated()
    {
        var user = new { Name = "John", Email = "john@example.com" };
        var content = new StringContent(
            JsonSerializer.Serialize(user),
            Encoding.UTF8,
            "application/json");

        var response = await _client.PostAsync("/api/users", content);

        response.StatusCode.Should().Be(HttpStatusCode.Created);
    }

    [Fact]
    public async Task CreateUser_InvalidUser_ReturnsBadRequest()
    {
        var user = new { Name = "" };  // Missing required fields
        var content = new StringContent(
            JsonSerializer.Serialize(user),
            Encoding.UTF8,
            "application/json");

        var response = await _client.PostAsync("/api/users", content);

        response.StatusCode.Should().Be(HttpStatusCode.BadRequest);
    }
}

// Custom factory with test configuration
public class CustomWebApplicationFactory : WebApplicationFactory<Program>
{
    protected override void ConfigureWebHost(IWebHostBuilder builder)
    {
        builder.ConfigureServices(services =>
        {
            // Replace database with in-memory
            var descriptor = services.SingleOrDefault(
                d => d.ServiceType == typeof(DbContextOptions<AppDbContext>));
            services.Remove(descriptor);

            services.AddDbContext<AppDbContext>(options =>
            {
                options.UseInMemoryDatabase("TestDb");
            });
        });
    }
}
```

---

## Configuration

### Test Project (.csproj)
```xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
    <IsPackable>false</IsPackable>
  </PropertyGroup>

  <ItemGroup>
    <!-- xUnit -->
    <PackageReference Include="xunit" Version="2.6.1" />
    <PackageReference Include="xunit.runner.visualstudio" Version="2.5.3" />
    
    <!-- Mocking -->
    <PackageReference Include="Moq" Version="4.20.69" />
    
    <!-- Assertions -->
    <PackageReference Include="FluentAssertions" Version="6.12.0" />
    
    <!-- Integration Tests -->
    <PackageReference Include="Microsoft.AspNetCore.Mvc.Testing" Version="8.0.0" />
    
    <!-- Coverage -->
    <PackageReference Include="coverlet.collector" Version="6.0.0" />
  </ItemGroup>

  <ItemGroup>
    <ProjectReference Include="..\MyApp\MyApp.csproj" />
  </ItemGroup>
</Project>
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
