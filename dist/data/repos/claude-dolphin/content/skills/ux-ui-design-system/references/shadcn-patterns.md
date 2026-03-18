# shadcn/ui Best Practices

Component patterns, organization, and integration guidelines for shadcn/ui.

## Philosophy

shadcn/ui is **not a component library** you install as a dependency. It's a collection of re-usable components that you **copy into your project**. You own the code and can customize freely.

### Key Principles

1. **Copy, don't install** - Components live in your codebase
2. **Customize freely** - You own the code, modify as needed
3. **Use Radix primitives** - Accessibility built-in
4. **Tailwind styling** - Utility-first CSS
5. **Design tokens** - CSS variables for theming

## Project Organization

### Directory Structure

```
src/
├── components/
│   ├── ui/           # shadcn components (copied)
│   │   ├── button.tsx
│   │   ├── card.tsx
│   │   ├── dialog.tsx
│   │   └── ...
│   ├── forms/        # Form-specific compositions
│   ├── layout/       # Layout components
│   └── [feature]/    # Feature-specific components
├── lib/
│   └── utils.ts      # cn() helper and utilities
└── styles/
    └── globals.css   # CSS variables, base styles
```

### Configuration Files

**components.json** - shadcn configuration:
```json
{
  "$schema": "https://ui.shadcn.com/schema.json",
  "style": "new-york",
  "rsc": true,
  "tsx": true,
  "tailwind": {
    "config": "tailwind.config.ts",
    "css": "src/styles/globals.css",
    "baseColor": "neutral",
    "cssVariables": true
  },
  "aliases": {
    "components": "@/components",
    "utils": "@/lib/utils",
    "ui": "@/components/ui",
    "lib": "@/lib",
    "hooks": "@/hooks"
  }
}
```

## Component Patterns

### Using the cn() Helper

Always use `cn()` for conditional class merging:

```tsx
import { cn } from "@/lib/utils"

function MyComponent({ className, variant }) {
  return (
    <div className={cn(
      "base-classes",
      variant === "primary" && "primary-classes",
      className // Allow parent override
    )}>
      {/* content */}
    </div>
  )
}
```

### Extending shadcn Components

**Option 1: Wrapper Component**
```tsx
// components/forms/submit-button.tsx
import { Button, ButtonProps } from "@/components/ui/button"
import { Loader2 } from "lucide-react"

interface SubmitButtonProps extends ButtonProps {
  isLoading?: boolean
}

export function SubmitButton({
  isLoading,
  children,
  disabled,
  ...props
}: SubmitButtonProps) {
  return (
    <Button disabled={disabled || isLoading} {...props}>
      {isLoading && <Loader2 className="mr-2 h-4 w-4 animate-spin" />}
      {children}
    </Button>
  )
}
```

**Option 2: Modify Source (for permanent changes)**
```tsx
// components/ui/button.tsx
// Add new variant directly to buttonVariants
const buttonVariants = cva(
  "...",
  {
    variants: {
      variant: {
        default: "...",
        destructive: "...",
        // Add custom variant
        brand: "bg-brand-500 text-white hover:bg-brand-600",
      },
    },
  }
)
```

### Composition Patterns

**Card Composition:**
```tsx
import {
  Card,
  CardContent,
  CardDescription,
  CardFooter,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"

function FeatureCard({ title, description, children, actions }) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>{title}</CardTitle>
        <CardDescription>{description}</CardDescription>
      </CardHeader>
      <CardContent>{children}</CardContent>
      {actions && <CardFooter>{actions}</CardFooter>}
    </Card>
  )
}
```

**Dialog Composition:**
```tsx
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogFooter,
  DialogHeader,
  DialogTitle,
  DialogTrigger,
} from "@/components/ui/dialog"

function ConfirmDialog({ trigger, title, description, onConfirm }) {
  const [open, setOpen] = useState(false)

  return (
    <Dialog open={open} onOpenChange={setOpen}>
      <DialogTrigger asChild>{trigger}</DialogTrigger>
      <DialogContent>
        <DialogHeader>
          <DialogTitle>{title}</DialogTitle>
          <DialogDescription>{description}</DialogDescription>
        </DialogHeader>
        <DialogFooter>
          <Button variant="outline" onClick={() => setOpen(false)}>
            Cancel
          </Button>
          <Button onClick={() => { onConfirm(); setOpen(false); }}>
            Confirm
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
```

## Form Integration

### React Hook Form + shadcn

```tsx
import { useForm } from "react-hook-form"
import { zodResolver } from "@hookform/resolvers/zod"
import * as z from "zod"
import {
  Form,
  FormControl,
  FormDescription,
  FormField,
  FormItem,
  FormLabel,
  FormMessage,
} from "@/components/ui/form"
import { Input } from "@/components/ui/input"
import { Button } from "@/components/ui/button"

const formSchema = z.object({
  email: z.string().email("Invalid email address"),
  name: z.string().min(2, "Name must be at least 2 characters"),
})

function MyForm() {
  const form = useForm<z.infer<typeof formSchema>>({
    resolver: zodResolver(formSchema),
    defaultValues: { email: "", name: "" },
  })

  function onSubmit(values: z.infer<typeof formSchema>) {
    console.log(values)
  }

  return (
    <Form {...form}>
      <form onSubmit={form.handleSubmit(onSubmit)} className="space-y-6">
        <FormField
          control={form.control}
          name="email"
          render={({ field }) => (
            <FormItem>
              <FormLabel>Email</FormLabel>
              <FormControl>
                <Input placeholder="email@example.com" {...field} />
              </FormControl>
              <FormDescription>Your email address</FormDescription>
              <FormMessage />
            </FormItem>
          )}
        />
        <Button type="submit">Submit</Button>
      </form>
    </Form>
  )
}
```

## Theming & Design Tokens

### CSS Variables Structure

```css
/* globals.css */
@layer base {
  :root {
    /* Background & Foreground */
    --background: 0 0% 100%;
    --foreground: 0 0% 3.9%;

    /* Primary */
    --primary: 0 0% 9%;
    --primary-foreground: 0 0% 98%;

    /* Secondary */
    --secondary: 0 0% 96.1%;
    --secondary-foreground: 0 0% 9%;

    /* Muted */
    --muted: 0 0% 96.1%;
    --muted-foreground: 0 0% 45.1%;

    /* Accent */
    --accent: 0 0% 96.1%;
    --accent-foreground: 0 0% 9%;

    /* Destructive */
    --destructive: 0 84.2% 60.2%;
    --destructive-foreground: 0 0% 98%;

    /* Border & Input */
    --border: 0 0% 89.8%;
    --input: 0 0% 89.8%;
    --ring: 0 0% 3.9%;

    /* Radius */
    --radius: 0.5rem;
  }

  .dark {
    --background: 0 0% 3.9%;
    --foreground: 0 0% 98%;
    /* ... dark mode values */
  }
}
```

### Adding Custom Colors

```css
:root {
  /* Brand colors */
  --brand: 220 90% 56%;
  --brand-foreground: 0 0% 100%;

  /* Status colors */
  --success: 142 76% 36%;
  --success-foreground: 0 0% 100%;
  --warning: 38 92% 50%;
  --warning-foreground: 0 0% 0%;
}
```

```ts
// tailwind.config.ts
export default {
  theme: {
    extend: {
      colors: {
        brand: "hsl(var(--brand))",
        "brand-foreground": "hsl(var(--brand-foreground))",
        success: "hsl(var(--success))",
        warning: "hsl(var(--warning))",
      },
    },
  },
}
```

## MCP Integration

### Available MCP Tools

When shadcn MCP is configured, use these tools:

**search_items_in_registries**
- Search for components by name or description
- Check multiple registries (default, custom)

**view_items_in_registries**
- Get full component source code
- Inspect dependencies

**get_item_examples_from_registries**
- Find demo code and usage patterns
- See component in context

**get_add_command_for_items**
- Get npx command to install components
- Include all dependencies

### MCP Workflow

1. **Before creating custom component:**
   - Search shadcn registry for existing solution
   - Check installed components in project
   - Only create custom if no suitable match

2. **When extending components:**
   - View existing source code
   - Find examples of similar patterns
   - Follow established conventions

3. **When troubleshooting:**
   - Get fresh component code from registry
   - Compare with local modifications
   - Identify divergence

## Common Mistakes

### ❌ Hardcoding Colors
```tsx
// Bad
<div className="bg-blue-500 text-white">
```

### ✅ Using Design Tokens
```tsx
// Good
<div className="bg-primary text-primary-foreground">
```

### ❌ Not Using cn()
```tsx
// Bad - classes override instead of merge
<Button className="mt-4 bg-red-500">
```

### ✅ Using cn() for Merging
```tsx
// Good - proper class merging
<Button className={cn("mt-4", variant === "danger" && "bg-red-500")}>
```

### ❌ Modifying node_modules
```tsx
// Bad - changes lost on npm install
// editing node_modules/@radix-ui/...
```

### ✅ Modifying Copied Source
```tsx
// Good - you own components/ui/*
// Modify components/ui/button.tsx directly
```

## Checklist

- [ ] Using cn() for all conditional classes
- [ ] Design tokens instead of hardcoded colors
- [ ] Extending via wrapper components or source modification
- [ ] Consistent composition patterns (Card*, Dialog*, etc.)
- [ ] Form components with React Hook Form + Zod
- [ ] Accessibility props maintained when extending
- [ ] Dark mode support via CSS variables
- [ ] Search MCP before creating custom components
