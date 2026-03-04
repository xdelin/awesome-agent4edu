#!/usr/bin/env node
/**
 * Reactome MCP Server
 * Production-ready Model Context Protocol server for Reactome pathway and systems biology data
 *
 * Copyright (c) 2025 Augmented Nature
 * Licensed under MIT License - see LICENSE file for details
 *
 * Developed by Augmented Nature - https://augmentednature.ai
 * Advancing AI for Scientific Discovery
 *
 * This server provides comprehensive access to Reactome's curated biological pathway data
 * through the Model Context Protocol, enabling AI systems to perform systems biology analysis.
 */
import { Server } from '@modelcontextprotocol/sdk/server/index.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import {
  CallToolRequestSchema,
  ErrorCode,
  ListResourcesRequestSchema,
  ListResourceTemplatesRequestSchema,
  ListToolsRequestSchema,
  McpError,
  ReadResourceRequestSchema,
} from '@modelcontextprotocol/sdk/types.js';
import axios, { AxiosInstance } from 'axios';

// Type guards and validation functions
const isValidSearchArgs = (args: any): args is { query: string; type?: string; size?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.query === 'string' &&
    args.query.length > 0 &&
    (args.type === undefined || ['pathway', 'reaction', 'protein', 'complex', 'disease'].includes(args.type)) &&
    (args.size === undefined || (typeof args.size === 'number' && args.size > 0 && args.size <= 100))
  );
};

const isValidIdArgs = (args: any): args is { id: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.id === 'string' &&
    args.id.length > 0
  );
};

const isValidGeneArgs = (args: any): args is { gene: string; species?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.gene === 'string' &&
    args.gene.length > 0 &&
    (args.species === undefined || typeof args.species === 'string')
  );
};

const isValidDiseaseArgs = (args: any): args is { disease: string; size?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.disease === 'string' &&
    args.disease.length > 0 &&
    (args.size === undefined || (typeof args.size === 'number' && args.size > 0 && args.size <= 100))
  );
};

const isValidInteractionArgs = (args: any): args is { pathwayId: string; interactionType?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.pathwayId === 'string' &&
    args.pathwayId.length > 0 &&
    (args.interactionType === undefined || ['protein-protein', 'regulatory', 'catalysis', 'all'].includes(args.interactionType))
  );
};

class ReactomeServer {
  private server: Server;
  private apiClient: AxiosInstance;

  constructor() {
    this.server = new Server(
      {
        name: 'reactome-server',
        version: '1.0.0',
      },
      {
        capabilities: {
          resources: {},
          tools: {},
        },
      }
    );

    // Initialize Reactome Content Service API client
    this.apiClient = axios.create({
      baseURL: 'https://reactome.org/ContentService',
      timeout: 30000,
      headers: {
        'User-Agent': 'Reactome-MCP-Server/1.0.0',
        'Content-Type': 'application/json',
      },
    });

    this.setupResourceHandlers();
    this.setupToolHandlers();

    // Error handling
    this.server.onerror = (error: Error) => console.error('[MCP Error]', error);
    process.on('SIGINT', async () => {
      await this.server.close();
      process.exit(0);
    });
  }

  private setupResourceHandlers() {
    this.server.setRequestHandler(
      ListResourceTemplatesRequestSchema,
      async () => ({
        resourceTemplates: [
          {
            uriTemplate: 'reactome://pathway/{id}',
            name: 'Reactome pathway information',
            mimeType: 'application/json',
            description: 'Complete pathway information for a Reactome pathway ID',
          },
          {
            uriTemplate: 'reactome://reaction/{id}',
            name: 'Reactome reaction information',
            mimeType: 'application/json',
            description: 'Detailed reaction information for a Reactome reaction ID',
          },
          {
            uriTemplate: 'reactome://protein/{id}',
            name: 'Reactome protein information',
            mimeType: 'application/json',
            description: 'Protein details and pathway associations',
          },
          {
            uriTemplate: 'reactome://disease/{id}',
            name: 'Reactome disease pathways',
            mimeType: 'application/json',
            description: 'Disease-associated pathways and mechanisms',
          },
          {
            uriTemplate: 'reactome://search/{query}',
            name: 'Reactome search results',
            mimeType: 'application/json',
            description: 'Search results across pathways, reactions, and entities',
          },
        ],
      })
    );

    this.server.setRequestHandler(
      ReadResourceRequestSchema,
      async (request: any) => {
        const uri = request.params.uri;

        // Handle pathway info requests
        const pathwayMatch = uri.match(/^reactome:\/\/pathway\/([A-Z0-9-]+)$/);
        if (pathwayMatch) {
          const pathwayId = pathwayMatch[1];
          try {
            const response = await this.apiClient.get(`/data/query/${pathwayId}`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch pathway ${pathwayId}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle reaction info requests
        const reactionMatch = uri.match(/^reactome:\/\/reaction\/([A-Z0-9-]+)$/);
        if (reactionMatch) {
          const reactionId = reactionMatch[1];
          try {
            const response = await this.apiClient.get(`/data/query/${reactionId}`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch reaction ${reactionId}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        throw new McpError(
          ErrorCode.InvalidRequest,
          `Invalid URI format: ${uri}`
        );
      }
    );
  }

  private setupToolHandlers() {
    this.server.setRequestHandler(ListToolsRequestSchema, async () => ({
      tools: [
        {
          name: 'search_pathways',
          description: 'Search for biological pathways by name, description, or keywords',
          inputSchema: {
            type: 'object',
            properties: {
              query: { type: 'string', description: 'Search query (pathway name, process, keywords)' },
              type: {
                type: 'string',
                enum: ['pathway', 'reaction', 'protein', 'complex', 'disease'],
                description: 'Type of entity to search for (default: pathway)'
              },
              size: { type: 'number', description: 'Number of results to return (1-100, default: 20)', minimum: 1, maximum: 100 },
            },
            required: ['query'],
          },
        },
        {
          name: 'get_pathway_details',
          description: 'Get comprehensive information about a specific pathway',
          inputSchema: {
            type: 'object',
            properties: {
              id: { type: 'string', description: 'Reactome pathway stable identifier (e.g., R-HSA-68886)' },
            },
            required: ['id'],
          },
        },
        {
          name: 'find_pathways_by_gene',
          description: 'Find all pathways containing a specific gene or protein',
          inputSchema: {
            type: 'object',
            properties: {
              gene: { type: 'string', description: 'Gene symbol or UniProt ID (e.g., BRCA1, P04637)' },
              species: { type: 'string', description: 'Species name or taxon ID (default: Homo sapiens)' },
            },
            required: ['gene'],
          },
        },
        {
          name: 'find_pathways_by_disease',
          description: 'Find disease-associated pathways and mechanisms',
          inputSchema: {
            type: 'object',
            properties: {
              disease: { type: 'string', description: 'Disease name or DOID identifier' },
              size: { type: 'number', description: 'Number of pathways to return (1-100, default: 25)', minimum: 1, maximum: 100 },
            },
            required: ['disease'],
          },
        },
        {
          name: 'get_pathway_hierarchy',
          description: 'Get hierarchical structure and parent/child relationships for a pathway',
          inputSchema: {
            type: 'object',
            properties: {
              id: { type: 'string', description: 'Reactome pathway stable identifier' },
            },
            required: ['id'],
          },
        },
        {
          name: 'get_pathway_participants',
          description: 'Get all molecules (proteins, genes, compounds) participating in a pathway',
          inputSchema: {
            type: 'object',
            properties: {
              id: { type: 'string', description: 'Reactome pathway stable identifier' },
            },
            required: ['id'],
          },
        },
        {
          name: 'get_pathway_reactions',
          description: 'Get all biochemical reactions within a pathway',
          inputSchema: {
            type: 'object',
            properties: {
              id: { type: 'string', description: 'Reactome pathway stable identifier' },
            },
            required: ['id'],
          },
        },
        {
          name: 'get_protein_interactions',
          description: 'Get protein-protein interactions within pathways',
          inputSchema: {
            type: 'object',
            properties: {
              pathwayId: { type: 'string', description: 'Reactome pathway stable identifier' },
              interactionType: {
                type: 'string',
                enum: ['protein-protein', 'regulatory', 'catalysis', 'all'],
                description: 'Type of interactions to retrieve (default: all)'
              },
            },
            required: ['pathwayId'],
          },
        },
      ],
    }));

    this.server.setRequestHandler(CallToolRequestSchema, async (request: any) => {
      const { name, arguments: args } = request.params;

      switch (name) {
        case 'search_pathways':
          return this.handleSearchPathways(args);
        case 'get_pathway_details':
          return this.handleGetPathwayDetails(args);
        case 'find_pathways_by_gene':
          return this.handleFindPathwaysByGene(args);
        case 'find_pathways_by_disease':
          return this.handleFindPathwaysByDisease(args);
        case 'get_pathway_hierarchy':
          return this.handleGetPathwayHierarchy(args);
        case 'get_pathway_participants':
          return this.handleGetPathwayParticipants(args);
        case 'get_pathway_reactions':
          return this.handleGetPathwayReactions(args);
        case 'get_protein_interactions':
          return this.handleGetProteinInteractions(args);
        default:
          throw new McpError(
            ErrorCode.MethodNotFound,
            `Unknown tool: ${name}`
          );
      }
    });
  }

  private async handleSearchPathways(args: any) {
    if (!isValidSearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid search arguments');
    }

    try {
      const params: any = {
        query: args.query,
        cluster: true,
      };

      if (args.type) {
        params.types = args.type;
      }

      const response = await this.apiClient.get('/search/query', { params });

      // Extract entries from all result groups
      let allEntries: any[] = [];
      if (response.data.results) {
        for (const group of response.data.results) {
          if (group.entries) {
            allEntries = allEntries.concat(group.entries);
          }
        }
      }

      // Filter by type if specified
      if (args.type) {
        const typeFilter = args.type.toLowerCase();
        allEntries = allEntries.filter((entry: any) =>
          entry.exactType?.toLowerCase().includes(typeFilter) ||
          entry.typeName?.toLowerCase().includes(typeFilter)
        );
      }

      // Limit results if specified
      if (args.size) {
        allEntries = allEntries.slice(0, args.size);
      }

      const formattedResults = {
        query: args.query,
        totalResults: response.data.numberOfMatches || 0,
        returnedResults: allEntries.length,
        results: allEntries.map((item: any) => ({
          id: item.stId || item.id,
          name: item.name?.replace(/<[^>]*>/g, '') || 'Unknown', // Remove HTML tags
          type: item.exactType || item.typeName || 'Unknown',
          species: Array.isArray(item.species) ? item.species[0] : item.species || 'Unknown',
          description: (item.summation?.substring(0, 200) || 'No description available') + '...',
          url: `https://reactome.org/content/detail/${item.stId || item.id}`
        }))
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(formattedResults, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error searching pathways: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async resolvePathwayId(identifier: string): Promise<string | null> {
    // If it's already a stable identifier, return it
    if (identifier.match(/^R-[A-Z]{3}-\d+$/)) {
      return identifier;
    }

    // Search for the pathway by name
    try {
      const searchResponse = await this.apiClient.get('/search/query', {
        params: {
          query: identifier,
          types: 'Pathway',
          cluster: true
        }
      });

      if (searchResponse.data.results &&
          searchResponse.data.results.length > 0 &&
          searchResponse.data.results[0].entries &&
          searchResponse.data.results[0].entries.length > 0) {
        const resolvedId = searchResponse.data.results[0].entries[0].stId;
        return resolvedId;
      }
    } catch (error) {
      // Silently handle pathway resolution errors
    }

    return null;
  }

  private async handleGetPathwayDetails(args: any) {
    if (!isValidIdArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Pathway ID is required');
    }

    try {
      // Resolve pathway ID if it's a name
      const pathwayId = await this.resolvePathwayId(args.id);
      if (!pathwayId) {
        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                error: `No pathway found for identifier: ${args.id}`,
                suggestion: 'Try using a Reactome stable identifier (e.g., R-HSA-1640170) or search for the pathway first'
              }, null, 2),
            },
          ],
          isError: true,
        };
      }

      // Get basic pathway information
      const basicInfo = await this.apiClient.get(`/data/query/${pathwayId}`);

      // Try to get additional pathway data using alternative endpoints
      let components = null;
      let participants = null;

      try {
        // Try the events endpoint
        const eventsResponse = await this.apiClient.get(`/data/events/${pathwayId}`);
        components = eventsResponse.data;
      } catch (e) {
        // Ignore component fetch errors
      }

      try {
        // Try the participants endpoint with different format
        const participantsResponse = await this.apiClient.get(`/data/participants/${pathwayId}`);
        participants = participantsResponse.data;
      } catch (e) {
        // Ignore participants fetch errors
      }

      const pathwayDetails = {
        id: pathwayId,
        originalQuery: args.id,
        basicInfo: basicInfo.data,
        components: components || 'Components data not available via API',
        participants: participants || 'Participants data not available via API',
        url: `https://reactome.org/content/detail/${pathwayId}`,
        diagramUrl: `https://reactome.org/PathwayBrowser/#/${pathwayId}`,
        note: 'Some detailed information may not be available through the current API endpoints'
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(pathwayDetails, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error getting pathway details: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleFindPathwaysByGene(args: any) {
    if (!isValidGeneArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid gene arguments');
    }

    try {
      // First search for the gene/protein entity
      const searchResponse = await this.apiClient.get('/search/query', {
        params: {
          query: args.gene,
          types: 'Protein',
          cluster: true
        }
      });

      // Extract protein entries from result groups
      let proteinEntries: any[] = [];
      if (searchResponse.data.results) {
        for (const group of searchResponse.data.results) {
          if (group.typeName === 'Protein' && group.entries) {
            proteinEntries = proteinEntries.concat(group.entries);
          }
        }
      }

      // Filter by species if specified
      if (args.species) {
        proteinEntries = proteinEntries.filter((entry: any) =>
          Array.isArray(entry.species) ?
            entry.species.some((s: string) => s.toLowerCase().includes(args.species!.toLowerCase())) :
            entry.species?.toLowerCase().includes(args.species!.toLowerCase())
        );
      }

      if (proteinEntries.length === 0) {
        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                gene: args.gene,
                species: args.species || 'Homo sapiens',
                message: 'No protein entity found for this gene',
                pathways: []
              }, null, 2),
            },
          ],
        };
      }

      // Get the first matching protein
      const protein = proteinEntries[0];

      // Find pathways containing this protein
      const pathwaysResponse = await this.apiClient.get(`/data/pathways/low/entity/${protein.stId}`);

      const result = {
        gene: args.gene,
        protein: {
          id: protein.stId,
          name: protein.name,
          species: protein.species?.[0]?.name
        },
        pathwayCount: pathwaysResponse.data?.length || 0,
        pathways: pathwaysResponse.data?.map((pathway: any) => ({
          id: pathway.stId,
          name: pathway.name,
          species: pathway.species?.[0]?.name,
          url: `https://reactome.org/content/detail/${pathway.stId}`
        })) || []
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error finding pathways by gene: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleFindPathwaysByDisease(args: any) {
    if (!isValidDiseaseArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid disease arguments');
    }

    try {
      // Search for disease-related pathways
      const searchResponse = await this.apiClient.get('/search/query', {
        params: {
          query: args.disease,
          types: 'Pathway',
          cluster: true
        }
      });

      // Extract pathway entries from result groups
      let pathwayEntries: any[] = [];
      if (searchResponse.data.results) {
        for (const group of searchResponse.data.results) {
          if (group.typeName === 'Pathway' && group.entries) {
            pathwayEntries = pathwayEntries.concat(group.entries);
          }
        }
      }

      // Limit results if specified
      if (args.size) {
        pathwayEntries = pathwayEntries.slice(0, args.size);
      }

      const diseasePathways = {
        disease: args.disease,
        pathwayCount: pathwayEntries.length,
        pathways: pathwayEntries.map((pathway: any) => ({
          id: pathway.stId || pathway.id,
          name: pathway.name?.replace(/<[^>]*>/g, '') || 'Unknown',
          type: pathway.exactType || pathway.typeName || 'Unknown',
          species: Array.isArray(pathway.species) ? pathway.species[0] : pathway.species || 'Unknown',
          description: (pathway.summation?.substring(0, 200) || 'No description available') + '...',
          url: `https://reactome.org/content/detail/${pathway.stId || pathway.id}`
        }))
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(diseasePathways, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error finding pathways by disease: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleGetPathwayHierarchy(args: any) {
    if (!isValidIdArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Pathway ID is required');
    }

    try {
      // Resolve pathway ID if it's a name
      const pathwayId = await this.resolvePathwayId(args.id);
      if (!pathwayId) {
        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                error: `No pathway found for identifier: ${args.id}`,
                suggestion: 'Try using a Reactome stable identifier (e.g., R-HSA-1640170) or search for the pathway first'
              }, null, 2),
            },
          ],
          isError: true,
        };
      }

      // Get basic pathway information first
      const pathwayInfo = await this.apiClient.get(`/data/query/${pathwayId}`);

      // Try alternative endpoints for hierarchy
      let ancestors = [];
      let children = [];

      try {
        // Try to get orthologous events (related pathways)
        const orthologousResponse = await this.apiClient.get(`/data/orthologous/${pathwayId}/pathways`);
        ancestors = orthologousResponse.data || [];
      } catch (e) {
        // Try to extract hierarchy info from basic pathway data
        if (pathwayInfo.data.hasEvent) {
          children = pathwayInfo.data.hasEvent.map((event: any) => ({
            id: event.stId || event.dbId,
            name: event.displayName || event.name,
            type: event.schemaClass || 'Event'
          }));
        }
      }

      const hierarchy = {
        pathwayId: pathwayId,
        originalQuery: args.id,
        basicInfo: {
          name: pathwayInfo.data.displayName || pathwayInfo.data.name,
          type: pathwayInfo.data.schemaClass,
          species: pathwayInfo.data.species?.[0]?.displayName
        },
        ancestors: ancestors.length > 0 ? ancestors.slice(0, 10).map((ancestor: any) => ({
          id: ancestor.stId || ancestor.dbId,
          name: ancestor.displayName || ancestor.name,
          type: ancestor.schemaClass || 'Pathway'
        })) : 'Ancestor information not available via current API',
        children: children.length > 0 ? children.slice(0, 10) : 'Child pathway information not available via current API',
        note: 'Hierarchy data may be limited due to API endpoint availability'
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(hierarchy, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error getting pathway hierarchy: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleGetPathwayParticipants(args: any) {
    if (!isValidIdArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Pathway ID is required');
    }

    try {
      // Resolve pathway ID if it's a name
      const pathwayId = await this.resolvePathwayId(args.id);
      if (!pathwayId) {
        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                error: `No pathway found for identifier: ${args.id}`,
                suggestion: 'Try using a Reactome stable identifier (e.g., R-HSA-1640170) or search for the pathway first'
              }, null, 2),
            },
          ],
          isError: true,
        };
      }

      // Try alternative approaches for getting participants
      let participants = [];

      try {
        // Try the working endpoint for participating molecules
        const response = await this.apiClient.get(`/data/pathway/${pathwayId}/participatingMolecules`);
        participants = response.data || [];
      } catch (error1) {
        try {
          // Alternative: get pathway details and extract participants from events
          const pathwayDetails = await this.apiClient.get(`/data/query/${pathwayId}`);
          if (pathwayDetails.data.hasEvent) {
            participants = pathwayDetails.data.hasEvent.map((event: any) => ({
              stId: event.stId || event.dbId,
              name: event.displayName || event.name,
              schemaClass: event.schemaClass,
              identifier: event.stId || event.dbId
            }));
          }
        } catch (error2) {
          // Final fallback: return pathway basic info
          const basicInfo = await this.apiClient.get(`/data/query/${pathwayId}`);
          return {
            content: [
              {
                type: 'text',
                text: JSON.stringify({
                  pathwayId: pathwayId,
                  originalQuery: args.id,
                  basicInfo: {
                    name: basicInfo.data.displayName || basicInfo.data.name,
                    type: basicInfo.data.schemaClass,
                    species: basicInfo.data.species?.[0]?.displayName
                  },
                  participantCount: 0,
                  participants: 'Participants data not available via current API endpoints',
                  note: 'Use search_pathways or get_pathway_details for alternative pathway information'
                }, null, 2),
              },
            ],
          };
        }
      }

      const result = {
        pathwayId: pathwayId,
        originalQuery: args.id,
        participantCount: participants.length,
        participants: participants.slice(0, 50).map((participant: any) => ({
          id: participant.stId,
          name: participant.name || participant.displayName,
          type: participant.schemaClass,
          species: participant.species?.[0]?.name || participant.species?.[0]?.displayName,
          identifier: participant.identifier,
          url: `https://reactome.org/content/detail/${participant.stId}`
        }))
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error getting pathway participants: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleGetPathwayReactions(args: any) {
    if (!isValidIdArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Pathway ID is required');
    }

    try {
      // Resolve pathway ID if it's a name
      const pathwayId = await this.resolvePathwayId(args.id);
      if (!pathwayId) {
        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                error: `No pathway found for identifier: ${args.id}`,
                suggestion: 'Try using a Reactome stable identifier (e.g., R-HSA-1640170) or search for the pathway first'
              }, null, 2),
            },
          ],
          isError: true,
        };
      }

      const response = await this.apiClient.get(`/data/pathway/${pathwayId}/containedEvents`);

      // Filter for reactions only
      const reactions = response.data?.filter((event: any) =>
        event.schemaClass === 'Reaction' || event.schemaClass === 'BlackBoxEvent'
      ) || [];

      const pathwayReactions = {
        pathwayId: pathwayId,
        originalQuery: args.id,
        reactionCount: reactions.length,
        reactions: reactions.map((reaction: any) => ({
          id: reaction.stId,
          name: reaction.name,
          type: reaction.schemaClass,
          reversible: reaction.reversible,
          species: reaction.species?.[0]?.name,
          url: `https://reactome.org/content/detail/${reaction.stId}`
        }))
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(pathwayReactions, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error getting pathway reactions: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleGetProteinInteractions(args: any) {
    if (!isValidInteractionArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid interaction arguments');
    }

    try {
      // Resolve pathway ID if it's a name
      const pathwayId = await this.resolvePathwayId(args.pathwayId);
      if (!pathwayId) {
        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                error: `No pathway found for identifier: ${args.pathwayId}`,
                suggestion: 'Try using a Reactome stable identifier (e.g., R-HSA-1640170) or search for the pathway first'
              }, null, 2),
            },
          ],
          isError: true,
        };
      }

      // Get basic pathway information first
      const pathwayInfo = await this.apiClient.get(`/data/query/${pathwayId}`);

      // Try multiple approaches to get interaction data
      let proteins = [];
      let reactions = [];

      try {
        // Try to get participating molecules
        const participantsResponse = await this.apiClient.get(`/data/pathway/${pathwayId}/participatingMolecules`);
        proteins = participantsResponse.data?.filter((p: any) =>
          p.schemaClass === 'EntityWithAccessionedSequence' || p.schemaClass === 'Protein'
        ) || [];
      } catch (e1) {
        // Alternative: search for proteins related to this pathway
        try {
          const searchResponse = await this.apiClient.get('/search/query', {
            params: {
              query: pathwayInfo.data.displayName || pathwayInfo.data.name,
              types: 'Protein',
              cluster: true
            }
          });

          if (searchResponse.data.results) {
            for (const group of searchResponse.data.results) {
              if (group.typeName === 'Protein' && group.entries) {
                proteins = group.entries.slice(0, 10); // Limit to 10 proteins
              }
            }
          }
        } catch (e2) {
          // Final fallback: extract from pathway hasEvent
          if (pathwayInfo.data.hasEvent) {
            proteins = pathwayInfo.data.hasEvent
              .filter((event: any) => event.schemaClass?.includes('Protein') || event.schemaClass?.includes('Entity'))
              .slice(0, 5);
          }
        }
      }

      try {
        // Try to get pathway reactions
        const reactionsResponse = await this.apiClient.get(`/data/pathway/${pathwayId}/containedEvents`);
        reactions = reactionsResponse.data?.filter((event: any) =>
          event.schemaClass === 'Reaction'
        ) || [];
      } catch (e) {
        // Extract reactions from pathway events
        if (pathwayInfo.data.hasEvent) {
          reactions = pathwayInfo.data.hasEvent
            .filter((event: any) => event.schemaClass === 'Reaction')
            .slice(0, 10);
        }
      }

      const interactions = {
        pathwayId: pathwayId,
        originalQuery: args.pathwayId,
        basicInfo: {
          name: pathwayInfo.data.displayName || pathwayInfo.data.name,
          type: pathwayInfo.data.schemaClass,
          species: pathwayInfo.data.species?.[0]?.displayName
        },
        proteinCount: proteins.length,
        reactionCount: reactions.length,
        proteins: proteins.slice(0, 20).map((protein: any) => ({
          id: protein.stId || protein.dbId,
          name: protein.name || protein.displayName,
          type: protein.schemaClass,
          identifier: protein.identifier
        })),
        potentialInteractions: reactions.slice(0, 15).map((reaction: any) => ({
          reactionId: reaction.stId || reaction.dbId,
          reactionName: reaction.name || reaction.displayName,
          type: reaction.schemaClass,
          reversible: reaction.reversible
        })),
        note: "Protein interactions inferred from pathway components and reactions. For detailed molecular interactions, consider using specialized protein interaction databases.",
        analysisNote: args.interactionType !== 'all' ? `Filtered for ${args.interactionType} interactions` : 'Showing all available interaction types'
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(interactions, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error getting protein interactions: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  async run() {
    const transport = new StdioServerTransport();
    await this.server.connect(transport);
    console.error('Reactome MCP server running on stdio');
  }
}

const server = new ReactomeServer();
server.run().catch(console.error);
