import { createExternalContextEntry } from '@/shared/mention/types';

describe('createExternalContextEntry', () => {
  it('creates entry with all fields', () => {
    const entry = createExternalContextEntry('/home/user/project', 'project', 'My Project');
    expect(entry).toEqual({
      contextRoot: '/home/user/project',
      folderName: 'project',
      displayName: 'My Project',
      displayNameLower: 'my project',
    });
  });

  it('lowercases displayName for displayNameLower', () => {
    const entry = createExternalContextEntry('/root', 'folder', 'CamelCaseName');
    expect(entry.displayNameLower).toBe('camelcasename');
  });

  it('handles empty strings', () => {
    const entry = createExternalContextEntry('', '', '');
    expect(entry).toEqual({
      contextRoot: '',
      folderName: '',
      displayName: '',
      displayNameLower: '',
    });
  });

  it('preserves original displayName casing', () => {
    const entry = createExternalContextEntry('/path', 'dir', 'UPPER');
    expect(entry.displayName).toBe('UPPER');
    expect(entry.displayNameLower).toBe('upper');
  });

  it('handles unicode characters', () => {
    const entry = createExternalContextEntry('/path', 'dir', 'Über Café');
    expect(entry.displayName).toBe('Über Café');
    expect(entry.displayNameLower).toBe('über café');
  });
});
