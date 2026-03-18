// Add type for detailed author info
export interface DetailedAuthorInfo {
  name: string;
  personal_name?: string;
  birth_date?: string;
  death_date?: string;
  bio?: string | { type: string; value: string }; // Bio can be string or object
  alternate_names?: string[];
  links?: { title: string; url: string; type: { key: string } }[];
  photos?: number[]; // Array of cover IDs
  source_records?: string[];
  wikipedia?: string;
  key: string;
  remote_ids?: {
    amazon?: string;
    librarything?: string;
    viaf?: string;
    goodreads?: string;
    storygraph?: string;
    wikidata?: string;
    isni?: string;
  };
  latest_revision?: number;
  revision: number;
  created?: { type: string; value: string };
  last_modified: { type: string; value: string };
}
