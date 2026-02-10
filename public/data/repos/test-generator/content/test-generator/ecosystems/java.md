# Java Testing Ecosystem

Comprehensive reference for testing Java projects.

## Detection

**Manifest files:** `pom.xml` (Maven), `build.gradle` / `build.gradle.kts` (Gradle)
**Test frameworks:** JUnit 5, JUnit 4, TestNG, Spock

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `junit-jupiter` in deps | JUnit 5 |
| `junit:junit:4.*` in deps | JUnit 4 |
| `org.testng` in deps | TestNG |
| `org.spockframework` in deps | Spock |
| `spring-boot-starter-test` | Spring Boot Test (JUnit 5) |
| `mockito-core` in deps | Mockito |

---

## File Structure

### Naming Conventions
```
*Test.java          # Standard suffix
*Tests.java         # Alternative
*IT.java            # Integration tests
*Spec.java          # Spock specifications
```

### Directory Patterns

**Maven/Gradle standard:**
```
src/
├── main/
│   └── java/
│       └── com/example/
│           ├── User.java
│           └── UserService.java
└── test/
    └── java/
        └── com/example/
            ├── UserTest.java
            └── UserServiceTest.java
```

**With integration tests:**
```
src/
├── main/java/
├── test/java/                    # Unit tests
│   └── com/example/
│       └── UserServiceTest.java
└── integrationTest/java/         # Integration tests
    └── com/example/
        └── UserServiceIT.java
```

---

## JUnit 5 Patterns

### Basic Structure
```java
import org.junit.jupiter.api.*;
import static org.junit.jupiter.api.Assertions.*;

class UserServiceTest {

    private UserService userService;

    @BeforeAll
    static void setUpClass() {
        // Run once before all tests
    }

    @AfterAll
    static void tearDownClass() {
        // Run once after all tests
    }

    @BeforeEach
    void setUp() {
        userService = new UserService();
    }

    @AfterEach
    void tearDown() {
        // Cleanup after each test
    }

    @Test
    @DisplayName("Should return user when valid ID provided")
    void shouldReturnUserWhenValidId() {
        // Arrange
        Long userId = 1L;

        // Act
        User result = userService.findById(userId);

        // Assert
        assertNotNull(result);
        assertEquals(userId, result.getId());
    }

    @Test
    @DisplayName("Should throw exception when user not found")
    void shouldThrowWhenUserNotFound() {
        // Arrange
        Long invalidId = 999L;

        // Act & Assert
        assertThrows(UserNotFoundException.class, () -> {
            userService.findById(invalidId);
        });
    }
}
```

### Assertions
```java
import static org.junit.jupiter.api.Assertions.*;

// Equality
assertEquals(expected, actual);
assertEquals(expected, actual, "failure message");
assertNotEquals(unexpected, actual);

// Boolean
assertTrue(condition);
assertFalse(condition);
assertNull(value);
assertNotNull(value);

// Same reference
assertSame(expected, actual);
assertNotSame(unexpected, actual);

// Arrays
assertArrayEquals(expectedArray, actualArray);

// Iterables
assertIterableEquals(expectedList, actualList);

// Exceptions
Exception exception = assertThrows(IllegalArgumentException.class, () -> {
    service.process(null);
});
assertEquals("Input cannot be null", exception.getMessage());

// Does not throw
assertDoesNotThrow(() -> service.process("valid"));

// Timeout
assertTimeout(Duration.ofSeconds(2), () -> {
    slowOperation();
});

// Grouped assertions (all run even if one fails)
assertAll("user properties",
    () -> assertEquals("John", user.getFirstName()),
    () -> assertEquals("Doe", user.getLastName()),
    () -> assertNotNull(user.getEmail())
);
```

### Parameterized Tests
```java
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.*;

class CalculatorTest {

    @ParameterizedTest
    @ValueSource(ints = {1, 2, 3, 4, 5})
    void shouldBePositive(int number) {
        assertTrue(number > 0);
    }

    @ParameterizedTest
    @NullAndEmptySource
    @ValueSource(strings = {"  ", "\t", "\n"})
    void shouldRejectBlankStrings(String input) {
        assertFalse(validator.isValid(input));
    }

    @ParameterizedTest
    @CsvSource({
        "1, 2, 3",
        "0, 0, 0",
        "-1, 1, 0",
        "100, 200, 300"
    })
    void shouldAddNumbers(int a, int b, int expected) {
        assertEquals(expected, calculator.add(a, b));
    }

    @ParameterizedTest
    @MethodSource("provideStringsForIsBlank")
    void shouldDetectBlankStrings(String input, boolean expected) {
        assertEquals(expected, StringUtils.isBlank(input));
    }

    static Stream<Arguments> provideStringsForIsBlank() {
        return Stream.of(
            Arguments.of(null, true),
            Arguments.of("", true),
            Arguments.of("  ", true),
            Arguments.of("not blank", false)
        );
    }

    @ParameterizedTest
    @EnumSource(Month.class)
    void shouldHandleAllMonths(Month month) {
        assertNotNull(calendar.getDaysIn(month));
    }

    @ParameterizedTest
    @EnumSource(value = Month.class, names = {"APRIL", "JUNE", "SEPTEMBER", "NOVEMBER"})
    void shouldHave30Days(Month month) {
        assertEquals(30, calendar.getDaysIn(month));
    }
}
```

### Nested Tests
```java
@DisplayName("UserService")
class UserServiceTest {

    @Nested
    @DisplayName("when creating user")
    class WhenCreatingUser {

        @Test
        @DisplayName("should save valid user")
        void shouldSaveValidUser() {
            // test
        }

        @Test
        @DisplayName("should reject duplicate email")
        void shouldRejectDuplicateEmail() {
            // test
        }
    }

    @Nested
    @DisplayName("when finding user")
    class WhenFindingUser {

        @Test
        @DisplayName("should return user by ID")
        void shouldReturnUserById() {
            // test
        }
    }
}
```

### Conditional Tests
```java
@Test
@EnabledOnOs(OS.LINUX)
void onlyOnLinux() {
}

@Test
@DisabledOnOs(OS.WINDOWS)
void notOnWindows() {
}

@Test
@EnabledOnJre(JRE.JAVA_17)
void onlyOnJava17() {
}

@Test
@EnabledIfEnvironmentVariable(named = "ENV", matches = "staging")
void onlyOnStaging() {
}

@Test
@Disabled("Temporarily disabled until bug #123 is fixed")
void temporarilyDisabled() {
}
```

---

## Mockito Patterns

### Basic Mocking
```java
import org.mockito.Mock;
import org.mockito.InjectMocks;
import org.mockito.junit.jupiter.MockitoExtension;
import static org.mockito.Mockito.*;

@ExtendWith(MockitoExtension.class)
class UserServiceTest {

    @Mock
    private UserRepository userRepository;

    @Mock
    private EmailService emailService;

    @InjectMocks
    private UserService userService;

    @Test
    void shouldFindUser() {
        // Arrange
        User expectedUser = new User(1L, "John");
        when(userRepository.findById(1L)).thenReturn(Optional.of(expectedUser));

        // Act
        User result = userService.findById(1L);

        // Assert
        assertEquals("John", result.getName());
        verify(userRepository).findById(1L);
    }

    @Test
    void shouldThrowWhenNotFound() {
        // Arrange
        when(userRepository.findById(anyLong())).thenReturn(Optional.empty());

        // Act & Assert
        assertThrows(UserNotFoundException.class, () -> {
            userService.findById(999L);
        });
    }
}
```

### Stubbing
```java
// Return value
when(mock.method()).thenReturn(value);
when(mock.method()).thenReturn(first, second, third);  // Consecutive

// Throw exception
when(mock.method()).thenThrow(new RuntimeException("error"));
when(mock.method()).thenThrow(RuntimeException.class);

// Call real method
when(mock.method()).thenCallRealMethod();

// Custom answer
when(mock.method(anyString())).thenAnswer(invocation -> {
    String arg = invocation.getArgument(0);
    return arg.toUpperCase();
});

// Void methods
doNothing().when(mock).voidMethod();
doThrow(new RuntimeException()).when(mock).voidMethod();
```

### Argument Matchers
```java
// Any
when(repo.findById(anyLong())).thenReturn(user);
when(repo.findByName(anyString())).thenReturn(user);
when(repo.save(any(User.class))).thenReturn(user);

// Specific conditions
when(repo.findByAge(eq(25))).thenReturn(users);
when(repo.findByAge(gt(18))).thenReturn(adults);
when(repo.findByName(startsWith("J"))).thenReturn(users);
when(repo.findByName(matches("\\w+@\\w+"))).thenReturn(users);

// Null handling
when(repo.findByName(isNull())).thenReturn(null);
when(repo.findByName(isNotNull())).thenReturn(user);

// Argument capture
ArgumentCaptor<User> captor = ArgumentCaptor.forClass(User.class);
verify(repo).save(captor.capture());
User savedUser = captor.getValue();
assertEquals("John", savedUser.getName());
```

### Verification
```java
// Basic verification
verify(mock).method();
verify(mock, times(2)).method();
verify(mock, never()).method();
verify(mock, atLeast(1)).method();
verify(mock, atMost(3)).method();

// Order verification
InOrder inOrder = inOrder(mock1, mock2);
inOrder.verify(mock1).firstMethod();
inOrder.verify(mock2).secondMethod();

// No more interactions
verifyNoMoreInteractions(mock);
verifyNoInteractions(mock);
```

---

## Spring Boot Test Patterns

### Unit Test (No Spring Context)
```java
@ExtendWith(MockitoExtension.class)
class UserServiceTest {

    @Mock
    private UserRepository userRepository;

    @InjectMocks
    private UserService userService;

    @Test
    void shouldFindUser() {
        // Pure unit test, no Spring
    }
}
```

### Integration Test (Full Context)
```java
@SpringBootTest
class UserServiceIntegrationTest {

    @Autowired
    private UserService userService;

    @Autowired
    private UserRepository userRepository;

    @Test
    void shouldCreateAndFindUser() {
        User user = userService.create("John", "john@example.com");
        User found = userService.findById(user.getId());
        assertEquals("John", found.getName());
    }
}
```

### Web Layer Test (MockMvc)
```java
@WebMvcTest(UserController.class)
class UserControllerTest {

    @Autowired
    private MockMvc mockMvc;

    @MockBean
    private UserService userService;

    @Test
    void shouldReturnUser() throws Exception {
        User user = new User(1L, "John");
        when(userService.findById(1L)).thenReturn(user);

        mockMvc.perform(get("/api/users/1"))
            .andExpect(status().isOk())
            .andExpect(jsonPath("$.name").value("John"));
    }

    @Test
    void shouldReturn404WhenNotFound() throws Exception {
        when(userService.findById(anyLong()))
            .thenThrow(new UserNotFoundException());

        mockMvc.perform(get("/api/users/999"))
            .andExpect(status().isNotFound());
    }

    @Test
    void shouldCreateUser() throws Exception {
        User user = new User(1L, "John");
        when(userService.create(any())).thenReturn(user);

        mockMvc.perform(post("/api/users")
                .contentType(MediaType.APPLICATION_JSON)
                .content("{\"name\":\"John\",\"email\":\"john@example.com\"}"))
            .andExpect(status().isCreated())
            .andExpect(jsonPath("$.id").value(1));
    }

    @Test
    void shouldValidateInput() throws Exception {
        mockMvc.perform(post("/api/users")
                .contentType(MediaType.APPLICATION_JSON)
                .content("{\"name\":\"\"}"))
            .andExpect(status().isBadRequest())
            .andExpect(jsonPath("$.errors.name").exists());
    }
}
```

### Data Layer Test
```java
@DataJpaTest
class UserRepositoryTest {

    @Autowired
    private UserRepository userRepository;

    @Autowired
    private TestEntityManager entityManager;

    @Test
    void shouldFindByEmail() {
        // Arrange
        User user = new User("John", "john@example.com");
        entityManager.persistAndFlush(user);

        // Act
        Optional<User> found = userRepository.findByEmail("john@example.com");

        // Assert
        assertTrue(found.isPresent());
        assertEquals("John", found.get().getName());
    }
}
```

### Test with Testcontainers
```java
@SpringBootTest
@Testcontainers
class UserServiceIntegrationTest {

    @Container
    static PostgreSQLContainer<?> postgres = new PostgreSQLContainer<>("postgres:15")
        .withDatabaseName("testdb")
        .withUsername("test")
        .withPassword("test");

    @DynamicPropertySource
    static void configureProperties(DynamicPropertyRegistry registry) {
        registry.add("spring.datasource.url", postgres::getJdbcUrl);
        registry.add("spring.datasource.username", postgres::getUsername);
        registry.add("spring.datasource.password", postgres::getPassword);
    }

    @Autowired
    private UserService userService;

    @Test
    void shouldPersistUser() {
        User user = userService.create("John", "john@example.com");
        assertNotNull(user.getId());
    }
}
```

---

## AssertJ (Fluent Assertions)

```java
import static org.assertj.core.api.Assertions.*;

// Strings
assertThat(name).isEqualTo("John");
assertThat(name).startsWith("Jo").endsWith("hn");
assertThat(name).containsIgnoringCase("OHN");

// Numbers
assertThat(age).isPositive();
assertThat(age).isBetween(18, 65);
assertThat(price).isCloseTo(10.0, within(0.01));

// Collections
assertThat(users).hasSize(3);
assertThat(users).contains(user1, user2);
assertThat(users).extracting("name").contains("John", "Jane");
assertThat(users).filteredOn("active", true).hasSize(2);

// Objects
assertThat(user)
    .hasFieldOrPropertyWithValue("name", "John")
    .hasFieldOrPropertyWithValue("email", "john@example.com");

// Exceptions
assertThatThrownBy(() -> service.process(null))
    .isInstanceOf(IllegalArgumentException.class)
    .hasMessageContaining("cannot be null");

assertThatCode(() -> service.process("valid"))
    .doesNotThrowAnyException();
```

---

## Configuration

### Maven (pom.xml)
```xml
<dependencies>
    <dependency>
        <groupId>org.junit.jupiter</groupId>
        <artifactId>junit-jupiter</artifactId>
        <version>5.10.0</version>
        <scope>test</scope>
    </dependency>
    <dependency>
        <groupId>org.mockito</groupId>
        <artifactId>mockito-junit-jupiter</artifactId>
        <version>5.5.0</version>
        <scope>test</scope>
    </dependency>
    <dependency>
        <groupId>org.assertj</groupId>
        <artifactId>assertj-core</artifactId>
        <version>3.24.2</version>
        <scope>test</scope>
    </dependency>
</dependencies>
```

### Gradle (build.gradle)
```groovy
dependencies {
    testImplementation 'org.junit.jupiter:junit-jupiter:5.10.0'
    testImplementation 'org.mockito:mockito-junit-jupiter:5.5.0'
    testImplementation 'org.assertj:assertj-core:3.24.2'
}

test {
    useJUnitPlatform()
}
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
