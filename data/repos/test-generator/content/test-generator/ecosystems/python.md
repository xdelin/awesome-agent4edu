# Python Testing Ecosystem

Comprehensive reference for testing Python projects.

## Detection

**Manifest files:** `pyproject.toml`, `requirements.txt`, `setup.py`, `setup.cfg`
**Test frameworks:** pytest (preferred), unittest, nose2

## Framework Detection

| Indicator | Framework |
|-----------|-----------|
| `pytest` in deps, `pytest.ini`, `pyproject.toml [tool.pytest]` | pytest |
| `conftest.py` present | pytest |
| Standard library usage, `unittest.TestCase` | unittest |
| `django.test` imports | Django Test |
| `flask.testing` imports | Flask Test |

---

## File Structure

### Naming Conventions
```
test_*.py          # pytest default (preferred)
*_test.py          # pytest alternative
tests.py           # Django default
```

### Directory Patterns

**Separate tests directory (recommended):**
```
project/
├── src/
│   └── mypackage/
│       ├── __init__.py
│       ├── services.py
│       └── utils.py
├── tests/
│   ├── __init__.py
│   ├── conftest.py          # Shared fixtures
│   ├── test_services.py
│   └── test_utils.py
```

**Tests alongside code:**
```
mypackage/
├── __init__.py
├── services.py
├── test_services.py
├── utils.py
└── test_utils.py
```

**Django pattern:**
```
myapp/
├── models.py
├── views.py
├── tests/
│   ├── __init__.py
│   ├── test_models.py
│   └── test_views.py
```

---

## pytest Patterns

### Basic Structure
```python
import pytest
from mypackage import function_to_test


class TestClassName:
    """Group related tests."""
    
    def setup_method(self):
        """Run before each test method."""
        self.resource = create_resource()
    
    def teardown_method(self):
        """Run after each test method."""
        cleanup_resource(self.resource)
    
    def test_should_return_expected_when_valid_input(self):
        """Descriptive test name."""
        # Arrange
        input_value = "test"
        
        # Act
        result = function_to_test(input_value)
        
        # Assert
        assert result == "expected"
    
    def test_should_raise_when_invalid_input(self):
        """Test exception handling."""
        with pytest.raises(ValueError, match="Invalid"):
            function_to_test(None)


# Alternative: function-based tests
def test_standalone_function():
    """Test without class."""
    assert function_to_test("input") == "expected"
```

### Assertions
```python
# Basic assertions
assert value == expected
assert value != unexpected
assert value is None
assert value is not None
assert value  # truthy
assert not value  # falsy

# Comparisons
assert value > 5
assert value >= 5
assert value < 10
assert value <= 10

# Membership
assert item in collection
assert item not in collection
assert "substring" in string

# Type checking
assert isinstance(value, str)
assert isinstance(value, (int, float))

# Approximate equality (floats)
assert value == pytest.approx(0.3, rel=1e-3)
assert value == pytest.approx(0.3, abs=0.01)

# Exception testing
with pytest.raises(ValueError):
    function_that_raises()

with pytest.raises(ValueError, match=r"must be positive"):
    function_that_raises()

# Warnings
with pytest.warns(DeprecationWarning):
    deprecated_function()
```

### Fixtures
```python
import pytest


@pytest.fixture
def sample_user():
    """Create a test user."""
    return User(id=1, name="Test User", email="test@example.com")


@pytest.fixture
def database_connection():
    """Provide database connection with cleanup."""
    conn = create_connection()
    yield conn
    conn.close()


@pytest.fixture(scope="module")
def expensive_resource():
    """Shared across all tests in module."""
    return create_expensive_resource()


@pytest.fixture(autouse=True)
def reset_state():
    """Automatically used by all tests."""
    global_state.reset()
    yield
    global_state.reset()


# Using fixtures
def test_user_name(sample_user):
    assert sample_user.name == "Test User"


def test_database_query(database_connection):
    result = database_connection.query("SELECT 1")
    assert result == 1
```

### conftest.py (Shared Fixtures)
```python
# tests/conftest.py
import pytest


@pytest.fixture(scope="session")
def app():
    """Create application instance for testing."""
    from myapp import create_app
    return create_app(testing=True)


@pytest.fixture
def client(app):
    """Test client for HTTP requests."""
    return app.test_client()


@pytest.fixture
def authenticated_client(client, sample_user):
    """Client with authentication."""
    client.post("/login", data={
        "email": sample_user.email,
        "password": "password"
    })
    return client
```

### Parametrized Tests
```python
import pytest


@pytest.mark.parametrize("input,expected", [
    ("hello", "HELLO"),
    ("world", "WORLD"),
    ("", ""),
    ("123", "123"),
])
def test_uppercase(input, expected):
    assert input.upper() == expected


@pytest.mark.parametrize("a,b,expected", [
    (1, 2, 3),
    (0, 0, 0),
    (-1, 1, 0),
    (100, 200, 300),
])
def test_add(a, b, expected):
    assert add(a, b) == expected


# Multiple parameter sets
@pytest.mark.parametrize("x", [1, 2, 3])
@pytest.mark.parametrize("y", [10, 20])
def test_multiply(x, y):
    assert multiply(x, y) == x * y  # 6 test combinations
```

### Markers
```python
import pytest


@pytest.mark.slow
def test_slow_operation():
    """Mark as slow, can be skipped with -m 'not slow'."""
    pass


@pytest.mark.skip(reason="Not implemented yet")
def test_future_feature():
    pass


@pytest.mark.skipif(sys.version_info < (3, 10), reason="Requires Python 3.10+")
def test_new_syntax():
    pass


@pytest.mark.xfail(reason="Known bug, see issue #123")
def test_known_failure():
    assert buggy_function() == "expected"
```

---

## Mocking

### unittest.mock
```python
from unittest.mock import Mock, patch, MagicMock


# Basic mock
mock = Mock()
mock.return_value = "result"
mock.method.return_value = "method result"

result = mock()
assert result == "result"
assert mock.called
assert mock.call_count == 1


# Patching
@patch("mypackage.module.external_service")
def test_with_patch(mock_service):
    mock_service.fetch.return_value = {"data": "mocked"}
    
    result = function_that_uses_service()
    
    assert result == {"data": "mocked"}
    mock_service.fetch.assert_called_once()


# Context manager
def test_with_context_patch():
    with patch("mypackage.module.external_service") as mock_service:
        mock_service.fetch.return_value = {"data": "mocked"}
        result = function_that_uses_service()
        assert result == {"data": "mocked"}


# Patch object attribute
@patch.object(MyClass, "method")
def test_patch_method(mock_method):
    mock_method.return_value = "mocked"
    instance = MyClass()
    assert instance.method() == "mocked"


# Side effects
mock.side_effect = ValueError("Error!")
mock.side_effect = [1, 2, 3]  # Return values in sequence
mock.side_effect = lambda x: x * 2  # Custom function


# Async mocking
from unittest.mock import AsyncMock

async_mock = AsyncMock(return_value="async result")
result = await async_mock()
```

### pytest-mock
```python
def test_with_mocker(mocker):
    """Using pytest-mock for simpler syntax."""
    mock_service = mocker.patch("mypackage.external_service")
    mock_service.fetch.return_value = {"data": "mocked"}
    
    result = function_that_uses_service()
    
    assert result == {"data": "mocked"}
```

---

## Django Testing

### Model Tests
```python
from django.test import TestCase
from myapp.models import User


class UserModelTest(TestCase):
    def setUp(self):
        self.user = User.objects.create(
            username="testuser",
            email="test@example.com"
        )
    
    def test_str_representation(self):
        self.assertEqual(str(self.user), "testuser")
    
    def test_email_unique(self):
        with self.assertRaises(IntegrityError):
            User.objects.create(
                username="another",
                email="test@example.com"  # Duplicate
            )
```

### View Tests
```python
from django.test import TestCase, Client
from django.urls import reverse


class ViewTests(TestCase):
    def setUp(self):
        self.client = Client()
        self.user = User.objects.create_user(
            username="testuser",
            password="testpass123"
        )
    
    def test_home_page(self):
        response = self.client.get(reverse("home"))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "home.html")
    
    def test_login_required(self):
        response = self.client.get(reverse("dashboard"))
        self.assertRedirects(response, "/login/?next=/dashboard/")
    
    def test_authenticated_access(self):
        self.client.login(username="testuser", password="testpass123")
        response = self.client.get(reverse("dashboard"))
        self.assertEqual(response.status_code, 200)
    
    def test_form_submission(self):
        self.client.login(username="testuser", password="testpass123")
        response = self.client.post(reverse("create"), {
            "title": "Test",
            "content": "Content"
        })
        self.assertEqual(response.status_code, 302)
```

### API Tests (Django REST Framework)
```python
from rest_framework.test import APITestCase
from rest_framework import status


class UserAPITests(APITestCase):
    def setUp(self):
        self.user = User.objects.create_user(
            username="testuser",
            password="testpass123"
        )
    
    def test_list_users(self):
        self.client.force_authenticate(user=self.user)
        response = self.client.get("/api/users/")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 1)
    
    def test_create_user(self):
        self.client.force_authenticate(user=self.user)
        response = self.client.post("/api/users/", {
            "username": "newuser",
            "email": "new@example.com"
        })
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
    
    def test_unauthorized(self):
        response = self.client.get("/api/users/")
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)
```

---

## Flask Testing

```python
import pytest
from myapp import create_app


@pytest.fixture
def app():
    app = create_app({"TESTING": True})
    yield app


@pytest.fixture
def client(app):
    return app.test_client()


def test_home(client):
    response = client.get("/")
    assert response.status_code == 200
    assert b"Welcome" in response.data


def test_json_api(client):
    response = client.post(
        "/api/users",
        json={"name": "Test", "email": "test@example.com"}
    )
    assert response.status_code == 201
    assert response.json["id"] is not None
```

---

## FastAPI Testing

```python
from fastapi.testclient import TestClient
from myapp import app


client = TestClient(app)


def test_read_root():
    response = client.get("/")
    assert response.status_code == 200
    assert response.json() == {"message": "Hello World"}


def test_create_item():
    response = client.post(
        "/items/",
        json={"name": "Test", "price": 10.5}
    )
    assert response.status_code == 201
    assert response.json()["name"] == "Test"


# Async tests
import pytest
from httpx import AsyncClient


@pytest.mark.anyio
async def test_async_endpoint():
    async with AsyncClient(app=app, base_url="http://test") as ac:
        response = await ac.get("/async-endpoint")
    assert response.status_code == 200
```

---

## Configuration

### pytest.ini
```ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = -v --tb=short
markers =
    slow: marks tests as slow
    integration: marks tests as integration tests
```

### pyproject.toml
```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
addopts = "-v --tb=short --strict-markers"

[tool.coverage.run]
source = ["src"]
omit = ["tests/*"]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "if TYPE_CHECKING:",
]
```

---

## Project-Specific Patterns

*This section grows with learnings specific to your project.*

---

## Learned Examples

*Real test examples from this project will be captured here.*
