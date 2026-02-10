# Rendrix - Master Your Commerce Universe

## Complete SaaS Development Plan & Implementation Guide

---

## Executive Summary

**Rendrix** is a comprehensive multi-tenant Software as a Service (SaaS) platform designed to empower entrepreneurs and businesses to create, manage, and scale multiple ecommerce stores from a single unified dashboard. Supporting diverse verticals—Toys, Kitchen Items, Nail Care, Home Decor, Garments, Beauty, Sports, Gadgets, and Home Appliances—Rendrix delivers enterprise-grade WooCommerce functionality with seamless multi-store orchestration.

---

## Development Prompt (Claude-Optimized)

```markdown
<role>
You are a senior full-stack architect and SaaS development expert specializing in
multi-tenant ecommerce platforms. You have deep expertise in:
- Multi-tenant architecture patterns (database-per-tenant, schema-per-tenant, shared-schema)
- WooCommerce ecosystem and WordPress headless implementations
- Payment gateway integrations (Stripe, PayPal, Square)
- Cloud infrastructure (AWS/GCP/Azure) with auto-scaling
- Modern frontend frameworks (Next.js, React) with TypeScript
- Database design for high-availability commerce systems
</role>

<context>
We are building "Rendrix - Master Your Commerce Universe", a SaaS platform enabling
users to create and manage multiple ecommerce websites across various business
verticals (Toy, Kitchen, Nail Care, Home, Garment, Beauty, Sports, Gadget,
Home Appliance). Each store requires full WooCommerce capabilities including
payment gateways, marketing tools, social media integration, customizable themes,
and SEO optimization. Users subscribe to plans (Free, Pro, Premium, Enterprise)
and manage all stores through a centralized dashboard with seamless store switching.
</context>

<task>
Design and implement the Rendrix SaaS platform following this phased approach:

Phase 1: Foundation Architecture
- Design multi-tenant database architecture with tenant isolation
- Implement authentication system with organization/workspace support
- Create subscription management with Stripe integration
- Build centralized dashboard shell with store context switching

Phase 2: Store Management Core
- Develop store provisioning and lifecycle management
- Implement WooCommerce headless API integration layer
- Create theme marketplace with preview and activation system
- Build SEO configuration module per store

Phase 3: Commerce Features
- Integrate payment gateways (Stripe, PayPal, Square, local options)
- Implement product catalog management with bulk operations
- Create order management and fulfillment workflows
- Build inventory synchronization across stores

Phase 4: Marketing & Growth
- Develop email marketing integration (Mailchimp, Klaviyo, SendGrid)
- Implement social media connector (Meta, TikTok, Pinterest)
- Create analytics dashboard with cross-store insights
- Build promotional tools (coupons, flash sales, bundles)

Phase 5: Enterprise Features
- Implement white-label capabilities
- Create API access for third-party integrations
- Build team collaboration with role-based access
- Develop custom domain management with SSL provisioning
</task>

<output_format>
Structure your implementation with:

<architecture>
Provide system architecture decisions with justifications
</architecture>

<database_schema>
Define complete database models with relationships
</database_schema>

<api_design>
Specify RESTful/GraphQL endpoints with request/response schemas
</api_design>

<implementation>
Deliver production-ready code with:
- TypeScript throughout
- Comprehensive error handling
- Input validation
- Security best practices
- Unit and integration tests
</implementation>

<deployment>
Include infrastructure-as-code and CI/CD configuration
</deployment>
</output_format>

<quality_requirements>
- Include as many relevant features as possible for a production-grade SaaS
- Apply security best practices: OWASP Top 10, encryption at rest/transit, audit logging
- Design for horizontal scalability supporting 10,000+ concurrent tenants
- Implement graceful degradation and circuit breaker patterns
- Follow 12-factor app methodology
- Ensure GDPR and PCI-DSS compliance readiness
</quality_requirements>

<examples>
<example name="store_switching">
// Seamless store context switching in the dashboard
const StoreSelector: React.FC = () => {
  const { stores, activeStore, switchStore } = useStoreContext();

  return (
    <Dropdown
      value={activeStore.id}
      options={stores.map(s => ({ value: s.id, label: s.name, icon: s.favicon }))}
      onChange={(storeId) => switchStore(storeId)}
      renderOption={(store) => (
        <StoreOption>
          <StoreFavicon src={store.icon} />
          <StoreName>{store.label}</StoreName>
          <StoreStatus status={store.status} />
        </StoreOption>
      )}
    />
  );
};
</example>

<example name="tenant_isolation">
// Middleware ensuring tenant data isolation
export const tenantIsolationMiddleware = async (req, res, next) => {
  const tenantId = req.user?.organizationId;
  if (!tenantId) return res.status(401).json({ error: 'Tenant context required' });

  // Attach tenant scope to all database queries
  req.db = createTenantScopedClient(tenantId);
  req.tenantId = tenantId;

  next();
};
</example>
</examples>
```

---

## System Architecture

### High-Level Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           RENDRIX ARCHITECTURE                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐                   │
│  │   Web App    │    │  Mobile App  │    │  Admin Panel │                   │
│  │  (Next.js)   │    │   (React     │    │   (Next.js)  │                   │
│  │              │    │    Native)   │    │              │                   │
│  └──────┬───────┘    └──────┬───────┘    └──────┬───────┘                   │
│         │                   │                   │                            │
│         └───────────────────┼───────────────────┘                            │
│                             │                                                │
│                    ┌────────▼────────┐                                       │
│                    │   API Gateway   │                                       │
│                    │  (Kong/AWS API) │                                       │
│                    └────────┬────────┘                                       │
│                             │                                                │
│         ┌───────────────────┼───────────────────┐                            │
│         │                   │                   │                            │
│  ┌──────▼──────┐    ┌───────▼───────┐   ┌──────▼──────┐                     │
│  │   Auth      │    │    Core API   │   │  Storefront │                     │
│  │  Service    │    │    Service    │   │    API      │                     │
│  │  (Auth0/    │    │   (Node.js)   │   │  (Node.js)  │                     │
│  │   Clerk)    │    │               │   │             │                     │
│  └──────┬──────┘    └───────┬───────┘   └──────┬──────┘                     │
│         │                   │                   │                            │
│         └───────────────────┼───────────────────┘                            │
│                             │                                                │
│              ┌──────────────┼──────────────┐                                 │
│              │              │              │                                 │
│       ┌──────▼──────┐ ┌─────▼─────┐ ┌──────▼──────┐                         │
│       │  PostgreSQL │ │   Redis   │ │ Elasticsearch│                        │
│       │  (Primary)  │ │  (Cache)  │ │  (Search)   │                         │
│       └─────────────┘ └───────────┘ └─────────────┘                         │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                        BACKGROUND SERVICES                           │    │
│  │  ┌───────────┐ ┌───────────┐ ┌───────────┐ ┌───────────┐            │    │
│  │  │  Queue    │ │  Email    │ │  Webhook  │ │   Cron    │            │    │
│  │  │  Worker   │ │  Service  │ │  Handler  │ │  Jobs     │            │    │
│  │  └───────────┘ └───────────┘ └───────────┘ └───────────┘            │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                      EXTERNAL INTEGRATIONS                           │    │
│  │  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐        │    │
│  │  │ Stripe  │ │ PayPal  │ │ Meta    │ │ Shippo  │ │ Mailchimp│       │    │
│  │  └─────────┘ └─────────┘ └─────────┘ └─────────┘ └─────────┘        │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Multi-Tenant Architecture Pattern

**Recommended: Hybrid Approach**

```
┌─────────────────────────────────────────────────────────────┐
│                    TENANT ISOLATION STRATEGY                 │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Shared Database + Schema-per-Tenant for Store Data         │
│                                                              │
│  ┌─────────────────────────────────────────────────────┐    │
│  │              SHARED TABLES (Platform)               │    │
│  │  • users, organizations, subscriptions              │    │
│  │  • plans, features, billing                         │    │
│  │  • themes_marketplace, integrations                 │    │
│  └─────────────────────────────────────────────────────┘    │
│                                                              │
│  ┌─────────────────────────────────────────────────────┐    │
│  │           TENANT-SCOPED TABLES (Store Data)         │    │
│  │  • stores, products, orders, customers              │    │
│  │  • inventory, coupons, reviews                      │    │
│  │  All queries filtered by tenant_id (RLS enforced)   │    │
│  └─────────────────────────────────────────────────────┘    │
│                                                              │
│  ┌─────────────────────────────────────────────────────┐    │
│  │           ISOLATED RESOURCES (Per Store)            │    │
│  │  • S3 buckets (media), CDN paths                    │    │
│  │  • Subdomain/custom domain routing                  │    │
│  │  • Separate Stripe Connect accounts                 │    │
│  └─────────────────────────────────────────────────────┘    │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

---

## Technology Stack

### Frontend

| Layer | Technology | Justification |
|-------|------------|---------------|
| Framework | Next.js 14 (App Router) | SSR/SSG, API routes, excellent DX |
| Language | TypeScript 5.x | Type safety, better tooling |
| Styling | Tailwind CSS + shadcn/ui | Rapid development, customizable |
| State | Zustand + TanStack Query | Lightweight, excellent caching |
| Forms | React Hook Form + Zod | Performant, schema validation |
| Charts | Recharts / Tremor | Dashboard analytics |

### Backend

| Layer | Technology | Justification |
|-------|------------|---------------|
| Runtime | Node.js 20 LTS | JavaScript ecosystem, async I/O |
| Framework | Fastify / NestJS | Performance, modularity |
| Language | TypeScript 5.x | Consistency with frontend |
| API | REST + GraphQL (hybrid) | Flexibility for different clients |
| Validation | Zod / class-validator | Runtime type checking |
| ORM | Prisma / Drizzle | Type-safe database access |

### Database & Storage

| Component | Technology | Justification |
|-----------|------------|---------------|
| Primary DB | PostgreSQL 16 | ACID, JSON support, RLS |
| Cache | Redis 7 | Session, queue, real-time |
| Search | Elasticsearch / Meilisearch | Product search, filtering |
| File Storage | AWS S3 / Cloudflare R2 | Media, exports, backups |
| CDN | Cloudflare / CloudFront | Global asset delivery |

### Infrastructure

| Component | Technology | Justification |
|-----------|------------|---------------|
| Cloud | AWS / GCP / Vercel | Scalability, managed services |
| Containers | Docker + Kubernetes | Orchestration, scaling |
| CI/CD | GitHub Actions | Automation, integration |
| Monitoring | Datadog / Grafana | Observability |
| Error Tracking | Sentry | Real-time error monitoring |

---

## Feature Modules

### 1. Authentication & Authorization

```typescript
// Feature: Multi-tenant Authentication System
interface AuthFeatures {
  authentication: {
    emailPassword: true;
    socialLogin: ['google', 'facebook', 'apple'];
    magicLink: true;
    twoFactorAuth: true;
  };

  authorization: {
    rbac: {
      roles: ['owner', 'admin', 'manager', 'staff', 'viewer'];
      permissions: string[]; // granular permissions
    };
    organizationScoping: true;
    storeScoping: true;
  };

  session: {
    jwtTokens: true;
    refreshTokenRotation: true;
    deviceManagement: true;
    concurrentSessionLimit: number;
  };
}
```

### 2. Subscription & Billing

```typescript
interface SubscriptionFeatures {
  plans: {
    free: {
      stores: 1;
      products: 50;
      bandwidth: '1GB';
      features: ['basic_themes', 'standard_support'];
    };
    pro: {
      stores: 3;
      products: 500;
      bandwidth: '10GB';
      features: ['premium_themes', 'seo_tools', 'priority_support'];
      price: { monthly: 29, yearly: 290 };
    };
    premium: {
      stores: 10;
      products: 5000;
      bandwidth: '100GB';
      features: ['all_themes', 'advanced_seo', 'marketing_suite', 'api_access'];
      price: { monthly: 79, yearly: 790 };
    };
    enterprise: {
      stores: 'unlimited';
      products: 'unlimited';
      bandwidth: 'unlimited';
      features: ['white_label', 'dedicated_support', 'custom_integrations', 'sla'];
      price: 'custom';
    };
  };

  billing: {
    provider: 'Stripe';
    features: ['invoicing', 'proration', 'dunning', 'usage_metering'];
  };
}
```

### 3. Store Management

```typescript
interface StoreManagement {
  provisioning: {
    instantSetup: true;
    templateStores: string[]; // industry-specific templates
    cloning: true;
    bulkImport: true;
  };

  configuration: {
    generalSettings: true;
    localization: {
      currencies: string[];
      languages: string[];
      timezones: string[];
    };
    taxConfiguration: true;
    shippingZones: true;
  };

  domains: {
    subdomain: true; // store-name.rendrix.com
    customDomain: true;
    sslProvisioning: 'automatic'; // Let's Encrypt
  };

  themes: {
    marketplace: true;
    customization: {
      colorSchemes: true;
      typography: true;
      layoutBuilder: true;
      customCSS: true;
      customJS: boolean; // premium only
    };
    preview: true;
    versionHistory: true;
  };
}
```

### 4. Product Catalog

```typescript
interface ProductCatalog {
  products: {
    types: ['simple', 'variable', 'grouped', 'digital', 'subscription'];
    attributes: true;
    variations: true;
    bulkEditing: true;
    import: ['csv', 'excel', 'woocommerce', 'shopify'];
    export: ['csv', 'excel'];
  };

  categories: {
    hierarchical: true;
    attributes: true;
    seoOptimized: true;
  };

  inventory: {
    tracking: true;
    lowStockAlerts: true;
    backorders: true;
    multiWarehouse: boolean; // enterprise
  };

  pricing: {
    regularPrice: true;
    salePrice: true;
    scheduledSales: true;
    tieredPricing: true;
    currencyConversion: true;
  };

  media: {
    images: { maxPerProduct: 20, optimization: true };
    videos: boolean;
    cdn: true;
    lazyLoading: true;
  };
}
```

### 5. Order Management

```typescript
interface OrderManagement {
  orders: {
    statuses: ['pending', 'processing', 'shipped', 'delivered', 'cancelled', 'refunded'];
    customStatuses: true;
    notes: true;
    timeline: true;
  };

  fulfillment: {
    pickingLists: true;
    packingSlips: true;
    shippingLabels: true;
    tracking: true;
    partialFulfillment: true;
  };

  returns: {
    rmaSystem: true;
    refundProcessing: true;
    restocking: true;
  };

  notifications: {
    email: true;
    sms: boolean;
    webhooks: true;
  };
}
```

### 6. Payment Integration

```typescript
interface PaymentIntegration {
  gateways: {
    stripe: { connect: true, direct: true };
    paypal: { standard: true, express: true };
    square: true;
    razorpay: true; // India
    mollie: true; // Europe
    afterpay: true; // BNPL
    klarna: true; // BNPL
  };

  features: {
    splitPayments: boolean; // marketplace
    recurringBilling: true;
    savedCards: true;
    multiCurrency: true;
    paymentLinks: true;
  };

  security: {
    pciCompliance: 'SAQ-A'; // tokenization
    fraudDetection: true;
    3dSecure: true;
  };
}
```

### 7. Marketing Suite

```typescript
interface MarketingSuite {
  email: {
    providers: ['mailchimp', 'klaviyo', 'sendgrid', 'built_in'];
    automation: {
      abandonedCart: true;
      welcomeSeries: true;
      winback: true;
      postPurchase: true;
    };
    templates: true;
    segmentation: true;
  };

  social: {
    metaCatalog: true; // Facebook/Instagram shops
    tiktokShop: true;
    pinterest: true;
    googleMerchant: true;
  };

  promotions: {
    coupons: {
      types: ['percentage', 'fixed', 'free_shipping', 'bogo'];
      conditions: true;
      usageLimits: true;
    };
    flashSales: true;
    bundles: true;
    loyaltyPoints: boolean;
  };

  analytics: {
    googleAnalytics: true;
    facebookPixel: true;
    customEvents: true;
    conversionTracking: true;
  };
}
```

### 8. SEO Module

```typescript
interface SEOModule {
  onPage: {
    metaTitles: true;
    metaDescriptions: true;
    urlStructure: 'customizable';
    canonicalUrls: true;
    schemaMarkup: {
      product: true;
      organization: true;
      breadcrumbs: true;
      reviews: true;
    };
  };

  technical: {
    sitemapGeneration: true;
    robotsTxt: true;
    redirectManager: true;
    pageSpeedOptimization: true;
  };

  tools: {
    keywordSuggestions: boolean;
    seoScoring: true;
    previewSnippets: true;
    bulkEditing: true;
  };
}
```

---

## Database Schema

### Core Entities

```sql
-- Organizations (Tenants)
CREATE TABLE organizations (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(100) UNIQUE NOT NULL,
    owner_id UUID NOT NULL REFERENCES users(id),
    subscription_id UUID REFERENCES subscriptions(id),
    settings JSONB DEFAULT '{}',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Users
CREATE TABLE users (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255),
    first_name VARCHAR(100),
    last_name VARCHAR(100),
    avatar_url TEXT,
    email_verified_at TIMESTAMPTZ,
    two_factor_enabled BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Organization Memberships
CREATE TABLE organization_members (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    organization_id UUID NOT NULL REFERENCES organizations(id) ON DELETE CASCADE,
    user_id UUID NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    role VARCHAR(50) NOT NULL DEFAULT 'member',
    permissions JSONB DEFAULT '[]',
    invited_at TIMESTAMPTZ,
    joined_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(organization_id, user_id)
);

-- Stores
CREATE TABLE stores (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    organization_id UUID NOT NULL REFERENCES organizations(id) ON DELETE CASCADE,
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(100) NOT NULL,
    subdomain VARCHAR(100) UNIQUE,
    custom_domain VARCHAR(255) UNIQUE,
    domain_verified BOOLEAN DEFAULT FALSE,
    industry VARCHAR(100),
    description TEXT,
    logo_url TEXT,
    favicon_url TEXT,
    theme_id UUID REFERENCES themes(id),
    theme_settings JSONB DEFAULT '{}',
    seo_settings JSONB DEFAULT '{}',
    settings JSONB DEFAULT '{}',
    status VARCHAR(50) DEFAULT 'active',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(organization_id, slug)
);

-- Enable Row Level Security
ALTER TABLE stores ENABLE ROW LEVEL SECURITY;

CREATE POLICY stores_tenant_isolation ON stores
    USING (organization_id = current_setting('app.current_tenant')::UUID);

-- Subscriptions
CREATE TABLE subscriptions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    organization_id UUID NOT NULL REFERENCES organizations(id) ON DELETE CASCADE,
    plan_id UUID NOT NULL REFERENCES plans(id),
    stripe_subscription_id VARCHAR(255),
    status VARCHAR(50) NOT NULL DEFAULT 'active',
    current_period_start TIMESTAMPTZ,
    current_period_end TIMESTAMPTZ,
    cancel_at TIMESTAMPTZ,
    canceled_at TIMESTAMPTZ,
    trial_ends_at TIMESTAMPTZ,
    metadata JSONB DEFAULT '{}',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Plans
CREATE TABLE plans (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name VARCHAR(100) NOT NULL,
    slug VARCHAR(50) UNIQUE NOT NULL,
    description TEXT,
    stripe_price_id_monthly VARCHAR(255),
    stripe_price_id_yearly VARCHAR(255),
    price_monthly DECIMAL(10,2),
    price_yearly DECIMAL(10,2),
    features JSONB NOT NULL DEFAULT '{}',
    limits JSONB NOT NULL DEFAULT '{}',
    is_active BOOLEAN DEFAULT TRUE,
    sort_order INTEGER DEFAULT 0,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- Products
CREATE TABLE products (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    store_id UUID NOT NULL REFERENCES stores(id) ON DELETE CASCADE,
    sku VARCHAR(100),
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(255) NOT NULL,
    description TEXT,
    short_description TEXT,
    type VARCHAR(50) DEFAULT 'simple',
    status VARCHAR(50) DEFAULT 'draft',
    visibility VARCHAR(50) DEFAULT 'visible',
    price DECIMAL(12,2),
    compare_at_price DECIMAL(12,2),
    cost_price DECIMAL(12,2),
    taxable BOOLEAN DEFAULT TRUE,
    tax_class VARCHAR(50),
    track_inventory BOOLEAN DEFAULT TRUE,
    quantity INTEGER DEFAULT 0,
    allow_backorders BOOLEAN DEFAULT FALSE,
    weight DECIMAL(10,3),
    dimensions JSONB,
    attributes JSONB DEFAULT '[]',
    metadata JSONB DEFAULT '{}',
    seo_title VARCHAR(255),
    seo_description TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    published_at TIMESTAMPTZ,
    UNIQUE(store_id, slug)
);

-- Product Variants
CREATE TABLE product_variants (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    product_id UUID NOT NULL REFERENCES products(id) ON DELETE CASCADE,
    sku VARCHAR(100),
    name VARCHAR(255),
    price DECIMAL(12,2),
    compare_at_price DECIMAL(12,2),
    quantity INTEGER DEFAULT 0,
    attributes JSONB NOT NULL DEFAULT '{}',
    image_url TEXT,
    sort_order INTEGER DEFAULT 0,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- Categories
CREATE TABLE categories (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    store_id UUID NOT NULL REFERENCES stores(id) ON DELETE CASCADE,
    parent_id UUID REFERENCES categories(id),
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(255) NOT NULL,
    description TEXT,
    image_url TEXT,
    sort_order INTEGER DEFAULT 0,
    seo_title VARCHAR(255),
    seo_description TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(store_id, slug)
);

-- Orders
CREATE TABLE orders (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    store_id UUID NOT NULL REFERENCES stores(id),
    order_number VARCHAR(50) NOT NULL,
    customer_id UUID REFERENCES customers(id),
    email VARCHAR(255) NOT NULL,
    phone VARCHAR(50),
    status VARCHAR(50) DEFAULT 'pending',
    payment_status VARCHAR(50) DEFAULT 'pending',
    fulfillment_status VARCHAR(50) DEFAULT 'unfulfilled',
    currency VARCHAR(3) DEFAULT 'USD',
    subtotal DECIMAL(12,2) NOT NULL,
    discount_total DECIMAL(12,2) DEFAULT 0,
    shipping_total DECIMAL(12,2) DEFAULT 0,
    tax_total DECIMAL(12,2) DEFAULT 0,
    total DECIMAL(12,2) NOT NULL,
    billing_address JSONB,
    shipping_address JSONB,
    notes TEXT,
    metadata JSONB DEFAULT '{}',
    placed_at TIMESTAMPTZ DEFAULT NOW(),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(store_id, order_number)
);

-- Order Items
CREATE TABLE order_items (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    order_id UUID NOT NULL REFERENCES orders(id) ON DELETE CASCADE,
    product_id UUID REFERENCES products(id),
    variant_id UUID REFERENCES product_variants(id),
    name VARCHAR(255) NOT NULL,
    sku VARCHAR(100),
    quantity INTEGER NOT NULL,
    price DECIMAL(12,2) NOT NULL,
    total DECIMAL(12,2) NOT NULL,
    metadata JSONB DEFAULT '{}'
);

-- Customers
CREATE TABLE customers (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    store_id UUID NOT NULL REFERENCES stores(id) ON DELETE CASCADE,
    email VARCHAR(255) NOT NULL,
    first_name VARCHAR(100),
    last_name VARCHAR(100),
    phone VARCHAR(50),
    accepts_marketing BOOLEAN DEFAULT FALSE,
    total_orders INTEGER DEFAULT 0,
    total_spent DECIMAL(12,2) DEFAULT 0,
    tags JSONB DEFAULT '[]',
    notes TEXT,
    metadata JSONB DEFAULT '{}',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(store_id, email)
);

-- Themes
CREATE TABLE themes (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(100) UNIQUE NOT NULL,
    description TEXT,
    preview_url TEXT,
    thumbnail_url TEXT,
    version VARCHAR(20),
    author VARCHAR(255),
    industries JSONB DEFAULT '[]',
    features JSONB DEFAULT '[]',
    is_premium BOOLEAN DEFAULT FALSE,
    price DECIMAL(10,2),
    settings_schema JSONB DEFAULT '{}',
    is_active BOOLEAN DEFAULT TRUE,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- Coupons
CREATE TABLE coupons (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    store_id UUID NOT NULL REFERENCES stores(id) ON DELETE CASCADE,
    code VARCHAR(50) NOT NULL,
    type VARCHAR(50) NOT NULL, -- percentage, fixed, free_shipping
    value DECIMAL(10,2),
    minimum_order DECIMAL(10,2),
    maximum_discount DECIMAL(10,2),
    usage_limit INTEGER,
    usage_count INTEGER DEFAULT 0,
    per_customer_limit INTEGER,
    applicable_products JSONB DEFAULT '[]',
    applicable_categories JSONB DEFAULT '[]',
    starts_at TIMESTAMPTZ,
    expires_at TIMESTAMPTZ,
    is_active BOOLEAN DEFAULT TRUE,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(store_id, code)
);

-- Indexes for Performance
CREATE INDEX idx_stores_organization ON stores(organization_id);
CREATE INDEX idx_products_store ON products(store_id);
CREATE INDEX idx_products_status ON products(store_id, status);
CREATE INDEX idx_orders_store ON orders(store_id);
CREATE INDEX idx_orders_status ON orders(store_id, status);
CREATE INDEX idx_orders_customer ON orders(customer_id);
CREATE INDEX idx_customers_store ON customers(store_id);
CREATE INDEX idx_customers_email ON customers(store_id, email);
```

---

## API Design

### RESTful Endpoints

```yaml
# Authentication
POST   /api/v1/auth/register
POST   /api/v1/auth/login
POST   /api/v1/auth/logout
POST   /api/v1/auth/refresh
POST   /api/v1/auth/forgot-password
POST   /api/v1/auth/reset-password
POST   /api/v1/auth/verify-email
POST   /api/v1/auth/2fa/enable
POST   /api/v1/auth/2fa/verify

# Organizations
GET    /api/v1/organizations
POST   /api/v1/organizations
GET    /api/v1/organizations/:id
PATCH  /api/v1/organizations/:id
DELETE /api/v1/organizations/:id
GET    /api/v1/organizations/:id/members
POST   /api/v1/organizations/:id/members/invite
DELETE /api/v1/organizations/:id/members/:userId

# Subscriptions
GET    /api/v1/subscriptions/plans
GET    /api/v1/subscriptions/current
POST   /api/v1/subscriptions
PATCH  /api/v1/subscriptions
DELETE /api/v1/subscriptions
POST   /api/v1/subscriptions/portal

# Stores (scoped to organization)
GET    /api/v1/stores
POST   /api/v1/stores
GET    /api/v1/stores/:storeId
PATCH  /api/v1/stores/:storeId
DELETE /api/v1/stores/:storeId
POST   /api/v1/stores/:storeId/clone

# Store Settings
GET    /api/v1/stores/:storeId/settings
PATCH  /api/v1/stores/:storeId/settings
GET    /api/v1/stores/:storeId/seo
PATCH  /api/v1/stores/:storeId/seo

# Themes
GET    /api/v1/themes
GET    /api/v1/themes/:themeId
GET    /api/v1/themes/:themeId/preview
POST   /api/v1/stores/:storeId/theme
GET    /api/v1/stores/:storeId/theme/settings
PATCH  /api/v1/stores/:storeId/theme/settings

# Products
GET    /api/v1/stores/:storeId/products
POST   /api/v1/stores/:storeId/products
GET    /api/v1/stores/:storeId/products/:productId
PATCH  /api/v1/stores/:storeId/products/:productId
DELETE /api/v1/stores/:storeId/products/:productId
POST   /api/v1/stores/:storeId/products/bulk
DELETE /api/v1/stores/:storeId/products/bulk
POST   /api/v1/stores/:storeId/products/import
GET    /api/v1/stores/:storeId/products/export

# Product Variants
GET    /api/v1/stores/:storeId/products/:productId/variants
POST   /api/v1/stores/:storeId/products/:productId/variants
PATCH  /api/v1/stores/:storeId/products/:productId/variants/:variantId
DELETE /api/v1/stores/:storeId/products/:productId/variants/:variantId

# Categories
GET    /api/v1/stores/:storeId/categories
POST   /api/v1/stores/:storeId/categories
GET    /api/v1/stores/:storeId/categories/:categoryId
PATCH  /api/v1/stores/:storeId/categories/:categoryId
DELETE /api/v1/stores/:storeId/categories/:categoryId

# Orders
GET    /api/v1/stores/:storeId/orders
POST   /api/v1/stores/:storeId/orders
GET    /api/v1/stores/:storeId/orders/:orderId
PATCH  /api/v1/stores/:storeId/orders/:orderId
POST   /api/v1/stores/:storeId/orders/:orderId/fulfill
POST   /api/v1/stores/:storeId/orders/:orderId/refund
POST   /api/v1/stores/:storeId/orders/:orderId/cancel

# Customers
GET    /api/v1/stores/:storeId/customers
POST   /api/v1/stores/:storeId/customers
GET    /api/v1/stores/:storeId/customers/:customerId
PATCH  /api/v1/stores/:storeId/customers/:customerId
DELETE /api/v1/stores/:storeId/customers/:customerId
GET    /api/v1/stores/:storeId/customers/:customerId/orders

# Coupons
GET    /api/v1/stores/:storeId/coupons
POST   /api/v1/stores/:storeId/coupons
GET    /api/v1/stores/:storeId/coupons/:couponId
PATCH  /api/v1/stores/:storeId/coupons/:couponId
DELETE /api/v1/stores/:storeId/coupons/:couponId

# Analytics
GET    /api/v1/stores/:storeId/analytics/overview
GET    /api/v1/stores/:storeId/analytics/sales
GET    /api/v1/stores/:storeId/analytics/products
GET    /api/v1/stores/:storeId/analytics/customers
GET    /api/v1/organizations/:id/analytics/cross-store

# Media
POST   /api/v1/stores/:storeId/media/upload
GET    /api/v1/stores/:storeId/media
DELETE /api/v1/stores/:storeId/media/:mediaId

# Storefront API (Public)
GET    /api/storefront/:storeSlug/products
GET    /api/storefront/:storeSlug/products/:slug
GET    /api/storefront/:storeSlug/categories
GET    /api/storefront/:storeSlug/categories/:slug
POST   /api/storefront/:storeSlug/cart
PATCH  /api/storefront/:storeSlug/cart
POST   /api/storefront/:storeSlug/checkout
POST   /api/storefront/:storeSlug/checkout/complete
```

---

## Development Phases

### Phase 1: Foundation (Weeks 1-4)

```
├── Infrastructure Setup
│   ├── Repository structure (monorepo with Turborepo)
│   ├── Docker development environment
│   ├── CI/CD pipelines (GitHub Actions)
│   ├── Database migrations setup (Prisma)
│   └── Environment configuration
│
├── Authentication System
│   ├── User registration/login
│   ├── JWT token management
│   ├── Password reset flow
│   ├── Email verification
│   └── OAuth providers (Google, GitHub)
│
├── Organization Management
│   ├── Organization CRUD
│   ├── Member invitations
│   ├── Role management
│   └── Organization switching
│
└── Subscription Foundation
    ├── Plan definitions
    ├── Stripe integration
    ├── Checkout flow
    └── Billing portal
```

### Phase 2: Store Management (Weeks 5-8)

```
├── Store Provisioning
│   ├── Store creation wizard
│   ├── Industry templates
│   ├── Subdomain assignment
│   └── Initial configuration
│
├── Theme System
│   ├── Theme marketplace UI
│   ├── Theme activation
│   ├── Settings customization
│   ├── Live preview
│   └── Color/typography controls
│
├── Domain Management
│   ├── Custom domain setup
│   ├── DNS verification
│   ├── SSL provisioning
│   └── Domain status monitoring
│
└── SEO Configuration
    ├── Meta tags management
    ├── URL structure settings
    ├── Sitemap generation
    └── Schema markup configuration
```

### Phase 3: Commerce Core (Weeks 9-14)

```
├── Product Management
│   ├── Product CRUD
│   ├── Variants and attributes
│   ├── Category management
│   ├── Bulk operations
│   ├── Import/export
│   └── Media management
│
├── Inventory System
│   ├── Stock tracking
│   ├── Low stock alerts
│   ├── Inventory history
│   └── Multi-location (Enterprise)
│
├── Order Management
│   ├── Order processing
│   ├── Status workflow
│   ├── Fulfillment
│   ├── Refunds/returns
│   └── Order notifications
│
└── Payment Integration
    ├── Stripe Connect
    ├── PayPal integration
    ├── Payment processing
    ├── Refund handling
    └── Payment webhooks
```

### Phase 4: Marketing & Growth (Weeks 15-18)

```
├── Email Marketing
│   ├── Provider integrations
│   ├── Automated flows
│   ├── Template builder
│   └── List management
│
├── Social Commerce
│   ├── Meta catalog sync
│   ├── Google Merchant
│   ├── Pinterest integration
│   └── Social proof widgets
│
├── Promotions
│   ├── Coupon system
│   ├── Flash sales
│   ├── Bundle builder
│   └── Discount rules engine
│
└── Analytics Dashboard
    ├── Sales analytics
    ├── Product performance
    ├── Customer insights
    └── Cross-store reporting
```

### Phase 5: Enterprise & Polish (Weeks 19-24)

```
├── Enterprise Features
│   ├── White-label branding
│   ├── API access
│   ├── Webhooks
│   ├── Advanced permissions
│   └── Audit logging
│
├── Performance Optimization
│   ├── Database optimization
│   ├── Caching strategy
│   ├── CDN configuration
│   └── Load testing
│
├── Security Hardening
│   ├── Security audit
│   ├── Penetration testing
│   ├── Compliance checks
│   └── Incident response
│
└── Launch Preparation
    ├── Documentation
    ├── Onboarding flows
    ├── Support systems
    └── Marketing site
```

---

## Security Considerations

### Authentication & Authorization

```typescript
// Security measures implementation
const securityConfig = {
  authentication: {
    passwordPolicy: {
      minLength: 12,
      requireUppercase: true,
      requireLowercase: true,
      requireNumbers: true,
      requireSymbols: true,
      preventCommon: true,
    },
    rateLimit: {
      login: { window: '15m', max: 5 },
      passwordReset: { window: '1h', max: 3 },
    },
    sessionManagement: {
      accessTokenExpiry: '15m',
      refreshTokenExpiry: '7d',
      rotateRefreshTokens: true,
      maxConcurrentSessions: 5,
    },
  },

  dataProtection: {
    encryptionAtRest: 'AES-256',
    encryptionInTransit: 'TLS 1.3',
    sensitiveFields: ['password', 'ssn', 'creditCard'],
    piiHandling: 'GDPR compliant',
  },

  apiSecurity: {
    cors: { strictOrigin: true },
    helmet: true,
    csrf: true,
    contentSecurityPolicy: true,
    rateLimiting: {
      global: { window: '1m', max: 100 },
      authenticated: { window: '1m', max: 500 },
    },
  },

  auditLogging: {
    events: [
      'auth.login',
      'auth.logout',
      'auth.passwordChange',
      'user.create',
      'user.delete',
      'store.create',
      'store.delete',
      'order.create',
      'payment.process',
      'settings.change',
    ],
    retention: '2 years',
  },
};
```

### Compliance Checklist

- [ ] **GDPR**: Data portability, right to deletion, consent management
- [ ] **PCI-DSS**: Tokenized payments, secure transmission, access controls
- [ ] **SOC 2**: Security policies, monitoring, incident response
- [ ] **CCPA**: California privacy rights compliance
- [ ] **OWASP Top 10**: Injection, XSS, CSRF, authentication, access control

---

## Deployment Architecture

```yaml
# docker-compose.production.yml
version: '3.8'

services:
  api:
    image: rendrix/api:latest
    deploy:
      replicas: 3
      resources:
        limits:
          cpus: '1'
          memory: 1G
    environment:
      - NODE_ENV=production
      - DATABASE_URL=${DATABASE_URL}
      - REDIS_URL=${REDIS_URL}
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:3000/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  web:
    image: rendrix/web:latest
    deploy:
      replicas: 2
    environment:
      - NEXT_PUBLIC_API_URL=${API_URL}

  worker:
    image: rendrix/worker:latest
    deploy:
      replicas: 2
    environment:
      - REDIS_URL=${REDIS_URL}

  postgres:
    image: postgres:16-alpine
    volumes:
      - postgres_data:/var/lib/postgresql/data
    environment:
      - POSTGRES_DB=rendrix
      - POSTGRES_USER=${DB_USER}
      - POSTGRES_PASSWORD=${DB_PASSWORD}

  redis:
    image: redis:7-alpine
    volumes:
      - redis_data:/data

volumes:
  postgres_data:
  redis_data:
```

---

## Success Metrics

| Metric | Target | Measurement |
|--------|--------|-------------|
| Platform Uptime | 99.9% | Monitoring (Datadog) |
| API Response Time (p95) | < 200ms | APM |
| Store Page Load Time | < 2s | Lighthouse/RUM |
| User Onboarding Completion | > 80% | Analytics |
| Monthly Churn Rate | < 3% | Billing data |
| NPS Score | > 50 | Surveys |
| Support Ticket Resolution | < 24h | Helpdesk |

---

## Next Steps

1. **Validate architecture** with technical review
2. **Set up development environment** with Docker
3. **Implement authentication** as foundation
4. **Build store provisioning** MVP
5. **Integrate Stripe** for subscriptions
6. **Launch beta** with early adopters

---

*Document Version: 1.0*
*Last Updated: December 2024*
*Prepared for: Rendrix Development Team*
