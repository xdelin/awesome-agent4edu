# PHP Testing Ecosystem

Comprehensive reference for testing PHP projects.

## Detection

**Manifest files:** `composer.json`
**Test frameworks:** PHPUnit, Pest, Codeception

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `phpunit/phpunit` in composer, `phpunit.xml` | PHPUnit |
| `pestphp/pest` in composer | Pest |
| `codeception/codeception` in composer | Codeception |
| `laravel/framework` in composer | Laravel (uses PHPUnit/Pest) |
| `symfony/framework-bundle` in composer | Symfony (uses PHPUnit) |

---

## File Structure

### Naming Conventions
```
*Test.php           # PHPUnit default
*.test.php          # Pest alternative
*Pest.php           # Pest explicit
```

### Directory Patterns

**Laravel pattern:**
```
tests/
├── Feature/                    # Integration/feature tests
│   ├── Auth/
│   │   ├── LoginTest.php
│   │   └── RegisterTest.php
│   └── Api/
│       └── UserControllerTest.php
├── Unit/                       # Unit tests
│   ├── Models/
│   │   └── UserTest.php
│   └── Services/
│       └── PaymentServiceTest.php
├── CreatesApplication.php
└── TestCase.php
```

**Standard PHP pattern:**
```
tests/
├── Unit/
│   └── CalculatorTest.php
├── Integration/
│   └── DatabaseTest.php
└── bootstrap.php
```

---

## PHPUnit Patterns

### Basic Structure
```php
<?php

namespace Tests\Unit;

use PHPUnit\Framework\TestCase;
use App\Calculator;

class CalculatorTest extends TestCase
{
    private Calculator $calculator;

    protected function setUp(): void
    {
        parent::setUp();
        $this->calculator = new Calculator();
    }

    protected function tearDown(): void
    {
        parent::tearDown();
    }

    public function test_should_add_two_numbers(): void
    {
        // Arrange
        $a = 5;
        $b = 3;

        // Act
        $result = $this->calculator->add($a, $b);

        // Assert
        $this->assertEquals(8, $result);
    }

    public function test_should_throw_on_division_by_zero(): void
    {
        $this->expectException(DivisionByZeroError::class);
        $this->expectExceptionMessage('Cannot divide by zero');

        $this->calculator->divide(10, 0);
    }
}
```

### Assertions
```php
// Equality
$this->assertEquals($expected, $actual);
$this->assertSame($expected, $actual);          // Strict type
$this->assertNotEquals($expected, $actual);

// Boolean
$this->assertTrue($value);
$this->assertFalse($value);
$this->assertNull($value);
$this->assertNotNull($value);

// Comparisons
$this->assertGreaterThan(5, $value);
$this->assertGreaterThanOrEqual(5, $value);
$this->assertLessThan(10, $value);

// Strings
$this->assertStringContainsString('needle', $haystack);
$this->assertStringStartsWith('prefix', $string);
$this->assertStringEndsWith('suffix', $string);
$this->assertMatchesRegularExpression('/pattern/', $string);

// Arrays
$this->assertCount(3, $array);
$this->assertContains($item, $array);
$this->assertArrayHasKey('key', $array);
$this->assertEmpty($array);
$this->assertNotEmpty($array);

// Objects
$this->assertInstanceOf(ClassName::class, $object);
$this->assertObjectHasProperty('property', $object);

// Files
$this->assertFileExists($path);
$this->assertFileEquals($expected, $actual);

// JSON
$this->assertJson($string);
$this->assertJsonStringEqualsJsonString($expected, $actual);

// Exceptions
$this->expectException(ExceptionClass::class);
$this->expectExceptionMessage('message');
$this->expectExceptionCode(404);
```

### Data Providers
```php
/**
 * @dataProvider additionProvider
 */
public function test_add(int $a, int $b, int $expected): void
{
    $this->assertEquals($expected, $this->calculator->add($a, $b));
}

public static function additionProvider(): array
{
    return [
        'positive numbers' => [1, 2, 3],
        'negative numbers' => [-1, -2, -3],
        'mixed numbers'    => [-1, 2, 1],
        'zeros'            => [0, 0, 0],
    ];
}


// Alternative: yield syntax
public static function validationProvider(): \Generator
{
    yield 'empty string' => ['', false];
    yield 'valid email' => ['test@example.com', true];
    yield 'invalid email' => ['not-an-email', false];
}
```

### Mocking
```php
use PHPUnit\Framework\MockObject\MockObject;

class ServiceTest extends TestCase
{
    private MockObject $repository;
    private UserService $service;

    protected function setUp(): void
    {
        $this->repository = $this->createMock(UserRepository::class);
        $this->service = new UserService($this->repository);
    }

    public function test_find_user(): void
    {
        $user = new User(1, 'John');

        $this->repository
            ->expects($this->once())
            ->method('find')
            ->with(1)
            ->willReturn($user);

        $result = $this->service->getUser(1);

        $this->assertEquals('John', $result->name);
    }

    public function test_throws_when_not_found(): void
    {
        $this->repository
            ->method('find')
            ->willReturn(null);

        $this->expectException(UserNotFoundException::class);

        $this->service->getUser(999);
    }

    public function test_consecutive_calls(): void
    {
        $this->repository
            ->method('find')
            ->willReturnOnConsecutiveCalls(
                new User(1, 'First'),
                new User(2, 'Second'),
                null
            );
    }
}
```

---

## Pest Patterns

### Basic Structure
```php
<?php

use App\Calculator;

beforeEach(function () {
    $this->calculator = new Calculator();
});

test('adds two numbers', function () {
    expect($this->calculator->add(2, 3))->toBe(5);
});

test('divides numbers', function () {
    expect($this->calculator->divide(10, 2))->toBe(5);
});

test('throws on division by zero', function () {
    $this->calculator->divide(10, 0);
})->throws(DivisionByZeroError::class);

it('can subtract numbers', function () {
    expect($this->calculator->subtract(5, 3))->toBe(2);
});
```

### Pest Expectations
```php
// Equality
expect($value)->toBe($expected);
expect($value)->toEqual($expected);
expect($value)->not->toBe($unexpected);

// Boolean
expect($value)->toBeTrue();
expect($value)->toBeFalse();
expect($value)->toBeNull();
expect($value)->toBeTruthy();
expect($value)->toBeFalsy();

// Types
expect($value)->toBeString();
expect($value)->toBeInt();
expect($value)->toBeFloat();
expect($value)->toBeArray();
expect($value)->toBeObject();
expect($value)->toBeInstanceOf(ClassName::class);

// Strings
expect($string)->toContain('substring');
expect($string)->toStartWith('prefix');
expect($string)->toEndWith('suffix');
expect($string)->toMatch('/pattern/');

// Arrays
expect($array)->toHaveCount(3);
expect($array)->toContain($item);
expect($array)->toHaveKey('key');
expect($array)->toMatchArray(['key' => 'value']);

// Chaining
expect($user)
    ->name->toBe('John')
    ->email->toContain('@')
    ->age->toBeGreaterThan(18);
```

### Pest Datasets
```php
dataset('emails', [
    'valid' => ['test@example.com', true],
    'invalid' => ['not-an-email', false],
    'empty' => ['', false],
]);

test('validates email', function (string $email, bool $expected) {
    expect(isValidEmail($email))->toBe($expected);
})->with('emails');

// Inline dataset
test('adds numbers', function (int $a, int $b, int $expected) {
    expect(add($a, $b))->toBe($expected);
})->with([
    [1, 2, 3],
    [0, 0, 0],
    [-1, 1, 0],
]);
```

---

## Laravel Testing

### Feature Tests (HTTP)
```php
<?php

namespace Tests\Feature;

use Tests\TestCase;
use App\Models\User;
use Illuminate\Foundation\Testing\RefreshDatabase;

class UserControllerTest extends TestCase
{
    use RefreshDatabase;

    public function test_can_list_users(): void
    {
        User::factory()->count(3)->create();

        $response = $this->getJson('/api/users');

        $response
            ->assertStatus(200)
            ->assertJsonCount(3, 'data');
    }

    public function test_can_create_user(): void
    {
        $response = $this->postJson('/api/users', [
            'name' => 'John Doe',
            'email' => 'john@example.com',
            'password' => 'password123',
        ]);

        $response
            ->assertStatus(201)
            ->assertJsonPath('data.name', 'John Doe');

        $this->assertDatabaseHas('users', [
            'email' => 'john@example.com',
        ]);
    }

    public function test_validates_required_fields(): void
    {
        $response = $this->postJson('/api/users', []);

        $response
            ->assertStatus(422)
            ->assertJsonValidationErrors(['name', 'email', 'password']);
    }

    public function test_requires_authentication(): void
    {
        $response = $this->getJson('/api/profile');

        $response->assertStatus(401);
    }

    public function test_authenticated_user_can_view_profile(): void
    {
        $user = User::factory()->create();

        $response = $this->actingAs($user)
            ->getJson('/api/profile');

        $response
            ->assertStatus(200)
            ->assertJsonPath('data.email', $user->email);
    }
}
```

### Validation Tests
```php
<?php

namespace Tests\Feature;

use Tests\TestCase;
use Illuminate\Foundation\Testing\RefreshDatabase;

class RegistrationValidationTest extends TestCase
{
    use RefreshDatabase;

    /**
     * @dataProvider invalidRegistrationData
     */
    public function test_validation_fails_with_invalid_data(
        array $data,
        array $expectedErrors
    ): void {
        $response = $this->postJson('/api/register', $data);

        $response
            ->assertStatus(422)
            ->assertJsonValidationErrors($expectedErrors);
    }

    public static function invalidRegistrationData(): array
    {
        return [
            'missing name' => [
                ['email' => 'test@example.com', 'password' => 'password123'],
                ['name'],
            ],
            'invalid email' => [
                ['name' => 'John', 'email' => 'not-an-email', 'password' => 'password123'],
                ['email'],
            ],
            'short password' => [
                ['name' => 'John', 'email' => 'test@example.com', 'password' => '123'],
                ['password'],
            ],
            'empty request' => [
                [],
                ['name', 'email', 'password'],
            ],
        ];
    }

    public function test_email_must_be_unique(): void
    {
        User::factory()->create(['email' => 'existing@example.com']);

        $response = $this->postJson('/api/register', [
            'name' => 'John',
            'email' => 'existing@example.com',
            'password' => 'password123',
        ]);

        $response->assertJsonValidationErrors(['email']);
    }
}
```

### Unit Tests (Models, Services)
```php
<?php

namespace Tests\Unit\Models;

use Tests\TestCase;
use App\Models\User;
use App\Models\Post;
use Illuminate\Foundation\Testing\RefreshDatabase;

class UserTest extends TestCase
{
    use RefreshDatabase;

    public function test_has_many_posts(): void
    {
        $user = User::factory()
            ->has(Post::factory()->count(3))
            ->create();

        $this->assertCount(3, $user->posts);
        $this->assertInstanceOf(Post::class, $user->posts->first());
    }

    public function test_full_name_attribute(): void
    {
        $user = User::factory()->create([
            'first_name' => 'John',
            'last_name' => 'Doe',
        ]);

        $this->assertEquals('John Doe', $user->full_name);
    }

    public function test_scope_active(): void
    {
        User::factory()->create(['active' => true]);
        User::factory()->create(['active' => false]);

        $activeUsers = User::active()->get();

        $this->assertCount(1, $activeUsers);
    }
}
```

### Mocking in Laravel
```php
<?php

namespace Tests\Feature;

use Tests\TestCase;
use App\Services\PaymentGateway;
use Mockery\MockInterface;

class PaymentTest extends TestCase
{
    public function test_processes_payment(): void
    {
        $this->mock(PaymentGateway::class, function (MockInterface $mock) {
            $mock->shouldReceive('charge')
                ->once()
                ->with(1000, 'tok_test')
                ->andReturn(['id' => 'ch_123', 'status' => 'succeeded']);
        });

        $response = $this->postJson('/api/payments', [
            'amount' => 1000,
            'token' => 'tok_test',
        ]);

        $response->assertStatus(200);
    }

    public function test_handles_payment_failure(): void
    {
        $this->mock(PaymentGateway::class, function (MockInterface $mock) {
            $mock->shouldReceive('charge')
                ->andThrow(new PaymentFailedException('Card declined'));
        });

        $response = $this->postJson('/api/payments', [
            'amount' => 1000,
            'token' => 'tok_test',
        ]);

        $response
            ->assertStatus(400)
            ->assertJsonPath('error', 'Card declined');
    }
}
```

### Browser Tests (Dusk)
```php
<?php

namespace Tests\Browser;

use Laravel\Dusk\Browser;
use Tests\DuskTestCase;
use App\Models\User;

class LoginTest extends DuskTestCase
{
    public function test_user_can_login(): void
    {
        $user = User::factory()->create();

        $this->browse(function (Browser $browser) use ($user) {
            $browser->visit('/login')
                ->type('email', $user->email)
                ->type('password', 'password')
                ->press('Login')
                ->assertPathIs('/dashboard')
                ->assertSee('Welcome');
        });
    }
}
```

---

## Configuration

### phpunit.xml
```xml
<?xml version="1.0" encoding="UTF-8"?>
<phpunit xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:noNamespaceSchemaLocation="vendor/phpunit/phpunit/phpunit.xsd"
         bootstrap="vendor/autoload.php"
         colors="true">
    <testsuites>
        <testsuite name="Unit">
            <directory>tests/Unit</directory>
        </testsuite>
        <testsuite name="Feature">
            <directory>tests/Feature</directory>
        </testsuite>
    </testsuites>
    <source>
        <include>
            <directory>app</directory>
        </include>
    </source>
    <php>
        <env name="APP_ENV" value="testing"/>
        <env name="DB_CONNECTION" value="sqlite"/>
        <env name="DB_DATABASE" value=":memory:"/>
    </php>
</phpunit>
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
