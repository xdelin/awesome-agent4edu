#!/usr/bin/env node

import { LinkedInMCPServer } from './server.js';
import { getConfig, validateConfig, needsOAuth } from './config.js';
import { OAuthManager } from './oauth-manager.js';
import { Logger } from './logger.js';
import type { TokenProvider } from './linkedin-client.js';

async function main() {
  try {
    const config = getConfig();
    validateConfig(config);

    let tokenProvider: TokenProvider | undefined;

    // If OAuth credentials are available, create an OAuthManager as the token
    // provider. It handles disk persistence (~/.config/linkedin-mcp/tokens.json),
    // automatic refresh via refresh_token, and falls back to the full OAuth
    // browser flow only when no valid token exists.
    if (needsOAuth(config)) {
      const logger = new Logger(config.logLevel);
      logger.info('Using OAuth token provider (tokens persisted to disk, auto-refresh enabled)');

      const oauthManager = new OAuthManager({
        clientId: config.linkedInClientId!,
        clientSecret: config.linkedInClientSecret!,
        redirectUri: config.linkedInRedirectUri!,
      }, logger);

      // Eagerly obtain a token so the OAuth browser flow (if needed) runs
      // before the MCP server starts accepting tool calls.
      await oauthManager.getAccessToken();

      // The OAuthManager implements TokenProvider — every API request will
      // call getAccessToken() which transparently refreshes expired tokens.
      tokenProvider = oauthManager;
    }

    const server = new LinkedInMCPServer(config, tokenProvider);
    await server.start();

    // Handle graceful shutdown
    process.on('SIGINT', async () => {
      await server.stop();
      process.exit(0);
    });

    process.on('SIGTERM', async () => {
      await server.stop();
      process.exit(0);
    });
  } catch (error) {
    console.error('Failed to start server:', error);
    process.exit(1);
  }
}

main();

