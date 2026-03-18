export interface OpenLibraryDoc {
  title: string;
  author_name?: string[];
  first_publish_year?: number;
  key: string; // Work key, e.g., "/works/OL45883W"
  edition_count: number;
  cover_i?: number; // Add optional cover ID
}

export interface OpenLibrarySearchResponse {
  docs: OpenLibraryDoc[];
}

export interface BookInfo {
  title: string;
  authors: string[];
  first_publish_year: number | null;
  open_library_work_key: string;
  edition_count: number;
  cover_url?: string;
}
