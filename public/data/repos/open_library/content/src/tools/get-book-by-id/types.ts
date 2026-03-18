interface OpenLibraryIdentifier {
  isbn_10?: string[];
  isbn_13?: string[];
  lccn?: string[];
  oclc?: string[];
  openlibrary?: string[];
}

interface OpenLibraryAuthor {
  url?: string;
  name: string;
}

interface OpenLibraryPublisher {
  name: string;
}

interface OpenLibraryCover {
  small?: string;
  medium?: string;
  large?: string;
}

// Represents the structure within the 'data' field of a record
interface OpenLibraryRecordData {
  url: string;
  key: string; // e.g., "/books/OL24194264M"
  title: string;
  subtitle?: string;
  authors?: OpenLibraryAuthor[];
  number_of_pages?: number;
  identifiers?: OpenLibraryIdentifier;
  publishers?: OpenLibraryPublisher[];
  publish_date?: string;
  subjects?: { name: string; url: string }[];
  ebooks?: { preview_url?: string; availability?: string; read_url?: string }[];
  cover?: OpenLibraryCover;
}

// Represents the structure within the 'details' field of a record
interface OpenLibraryRecordDetails {
  bib_key: string;
  info_url: string;
  preview?: string;
  preview_url?: string;
  thumbnail_url?: string;
  details?: {
    // Yes, there's another nested 'details' sometimes
    key: string;
    title: string;
    authors?: OpenLibraryAuthor[];
    publishers?: OpenLibraryPublisher[];
    publish_date?: string;
    works?: { key: string }[]; // e.g., "/works/OL15610910W"
    covers?: number[]; // Cover IDs, not URLs
    lccn?: string[];
    oclc_numbers?: string[];
    isbn_10?: string[];
    isbn_13?: string[];
    number_of_pages?: number;
    // ... other potential fields
  };
}

// Represents a single record in the 'records' object
export interface OpenLibraryRecord {
  recordURL: string;
  data: OpenLibraryRecordData;
  details: OpenLibraryRecordDetails; // This is the details object we need
  // ... other potential fields like isbns, olids etc. at this level
}

// Type for the overall API response structure
export interface OpenLibraryBookResponse {
  // The keys are dynamic (e.g., "/books/OL24194264M")
  records: Record<string, OpenLibraryRecord>;
  items: unknown[]; // Items might not be needed for core details
}

// --- Formatted Book Details returned by the tool --- //

export interface BookDetails {
  title: string;
  subtitle?: string;
  authors: string[];
  publishers?: string[];
  publish_date?: string;
  number_of_pages?: number;
  isbn_13?: string[];
  isbn_10?: string[];
  lccn?: string[];
  oclc?: string[];
  olid?: string[]; // Add OLID field
  open_library_edition_key: string; // e.g., "/books/OL24194264M"
  open_library_work_key?: string; // e.g., "/works/OL15610910W"
  cover_url?: string;
  info_url: string;
  preview_url?: string;
}
