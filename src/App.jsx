
import { useState, useEffect, useMemo } from 'react';
import skillsData from './data/skills.json';
import repoMap from './data/repo_map.json';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { oneLight } from 'react-syntax-highlighter/dist/esm/styles/prism';
import { 
  Search, ExternalLink, BookOpen, Code, PenTool, Database, Layout, 
  Brain, Calculator, Briefcase, Box, Menu, X, Filter, Globe, Loader2, 
  Info, Folder, File as FileIcon, ChevronRight, ChevronDown, FileCode, 
  FileJson, FileText, Image as ImageIcon
} from 'lucide-react';

const CATEGORY_ICONS = {
  // English mappings
  "Intelligent Tutoring": BookOpen,
  "Math & Science": Calculator,
  "Computer Science": Code,
  "Data & Analysis": Database,
  "Visual & Presentation": Layout,
  "Academic & Writing": PenTool,
  "Notes & Knowledge Base": Brain,
  "Career & Productivity": Briefcase,
  "Meta Skills": Box,
  // Chinese mappings
  "智能导学": BookOpen,
  "数理科学": Calculator,
  "计算机科学": Code,
  "数据与分析": Database,
  "视觉与演示 (PPT)": Layout,
  "学术与写作": PenTool,
  "笔记与知识库": Brain,
  "职业规划与生产力": Briefcase,
  "元技能": Box
};

// File Icon Helper
const getFileIcon = (filename) => {
  const ext = filename.split('.').pop().toLowerCase();
  switch(ext) {
    case 'js':
    case 'jsx':
    case 'ts':
    case 'tsx':
    case 'py':
    case 'html':
    case 'css':
      return <FileCode size={16} className="text-blue-500" />;
    case 'json':
    case 'yml':
    case 'yaml':
      return <FileJson size={16} className="text-yellow-500" />;
    case 'md':
    case 'txt':
      return <FileText size={16} className="text-slate-500" />;
    case 'png':
    case 'jpg':
    case 'jpeg':
    case 'svg':
      return <ImageIcon size={16} className="text-purple-500" />;
    default:
      return <FileIcon size={16} className="text-slate-400" />;
  }
};

const FileTreeNode = ({ node, level, onSelect, selectedPath }) => {
  const [isOpen, setIsOpen] = useState(false);
  const isSelected = selectedPath === node.path;

  // Auto-expand if selected file is inside this folder (simple heuristic: path starts with)
  useEffect(() => {
    if (selectedPath && selectedPath.startsWith(node.path + '/')) {
      setIsOpen(true);
    }
  }, [selectedPath, node.path]);

  if (node.type === 'folder') {
    return (
      <div>
        <div 
          className="flex items-center gap-1.5 py-1 px-2 hover:bg-stone-100 cursor-pointer rounded-md text-sm text-stone-700 select-none whitespace-nowrap"
          style={{ paddingLeft: `${level * 12 + 8}px` }}
          onClick={() => setIsOpen(!isOpen)}
        >
          {isOpen ? <ChevronDown size={14} className="text-stone-400" /> : <ChevronRight size={14} className="text-stone-400" />}
          <Folder size={16} className="text-orange-400 fill-orange-100" />
          <span className="truncate">{node.name}</span>
        </div>
        {isOpen && (
          <div>
            {node.children.map((child, idx) => (
              <FileTreeNode 
                key={idx} 
                node={child} 
                level={level + 1} 
                onSelect={onSelect}
                selectedPath={selectedPath}
              />
            ))}
          </div>
        )}
      </div>
    );
  }

  return (
    <div 
      className={`flex items-center gap-1.5 py-1 px-2 cursor-pointer rounded-md text-sm select-none whitespace-nowrap
        ${isSelected ? 'bg-orange-100 text-orange-900 font-medium' : 'text-stone-600 hover:bg-stone-100'}`}
      style={{ paddingLeft: `${level * 12 + 24}px` }}
      onClick={() => onSelect(node)}
    >
      {getFileIcon(node.name)}
      <span className="truncate">{node.name}</span>
    </div>
  );
};

function RepoBrowser({ skill, repoId, onClose }) {
  const [manifest, setManifest] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [selectedFile, setSelectedFile] = useState(null);
  const [fileContent, setFileContent] = useState(null);
  const [contentLoading, setContentLoading] = useState(false);

  // Initial load: Fetch manifest
  useEffect(() => {
    const fetchManifest = async () => {
      try {
        setLoading(true);
        const res = await fetch(`data/repos/${repoId}/manifest.json`);
        if (!res.ok) throw new Error('Repository data not available locally');
        const data = await res.json();
        setManifest(data);
        
        // Try to select README by default
        const findReadme = (nodes) => {
           for (let node of nodes) {
              if (node.type === 'file' && node.name.toLowerCase().startsWith('readme')) return node;
              if (node.type === 'folder') {
                 const found = findReadme(node.children);
                 if (found) return found;
              }
           }
           return null;
        };
        const readme = findReadme(data.tree);
        if (readme) handleFileSelect(readme);
        
      } catch (err) {
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };
    fetchManifest();
  }, [repoId]);

  const handleFileSelect = async (fileNode) => {
     setSelectedFile(fileNode);
     setContentLoading(true);
     setFileContent(null);
     
     try {
       // Check for image
       const isImage = /\.(png|jpg|jpeg|gif|svg|ico)$/i.test(fileNode.name);
       
       const res = await fetch(`data/repos/${repoId}/content/${fileNode.path}`);
       if (!res.ok) throw new Error("Failed to load file content");
       
       if (isImage) {
           const blob = await res.blob();
           setFileContent(URL.createObjectURL(blob));
       } else {
           const text = await res.text();
           setFileContent(text);
       }
     } catch (err) {
       setFileContent(`Error loading file: ${err.message}`);
     } finally {
       setContentLoading(false);
     }
  };

  const isImage = selectedFile && /\.(png|jpg|jpeg|gif|svg|ico)$/i.test(selectedFile.name);

  return (
    <div className="fixed inset-0 z-[60] flex items-center justify-center p-2 sm:p-4" role="dialog">
       <div className="fixed inset-0 bg-stone-900/60 backdrop-blur-sm" onClick={onClose} />
       
       <div className="relative bg-white rounded-xl shadow-2xl w-full max-w-[95vw] h-[90vh] flex flex-col overflow-hidden animate-in zoom-in-95 duration-200">
          
          {/* Header */}
          <div className="flex items-center justify-between px-4 py-3 border-b border-stone-200 bg-stone-50">
             <div className="flex items-center gap-3 overflow-hidden">
                <div className="p-1.5 bg-orange-100 rounded-lg text-orange-600 shrink-0">
                   <Code size={20} />
                </div>
                <div className="min-w-0">
                   <h3 className="text-base font-bold text-stone-900 truncate">{skill.name}</h3>
                   <p className="text-xs text-stone-500 truncate font-mono">/{repoId}</p>
                </div>
             </div>
             <div className="flex gap-2 shrink-0">
                <a href={skill.url} target="_blank" rel="noopener noreferrer" className="p-2 text-stone-500 hover:text-stone-900 hover:bg-stone-200 rounded-lg transition-colors">
                   <ExternalLink size={18} />
                </a>
                <button onClick={onClose} className="p-2 text-stone-500 hover:text-red-600 hover:bg-red-50 rounded-lg transition-colors">
                   <X size={18} />
                </button>
             </div>
          </div>

          {/* Main Body */}
          <div className="flex-1 flex min-h-0">
             
             {/* Sidebar: File Tree */}
             <div className="w-64 md:w-80 border-r border-stone-200 bg-stone-50 flex flex-col min-h-0">
                <div className="p-3 border-b border-stone-100">
                   <h4 className="text-xs font-bold text-stone-400 uppercase tracking-wider">Explorer</h4>
                </div>
                <div className="flex-1 overflow-y-auto custom-scrollbar p-2">
                   {loading && (
                      <div className="flex flex-col items-center justify-center py-10 text-stone-400 gap-2">
                         <Loader2 size={20} className="animate-spin" />
                         <span className="text-xs">Loading tree...</span>
                      </div>
                   )}
                   {error && <div className="text-xs text-red-500 p-4 text-center">{error}</div>}
                   
                   {manifest && manifest.tree.map((node, idx) => (
                      <FileTreeNode 
                         key={idx} 
                         node={node} 
                         level={0} 
                         onSelect={handleFileSelect} 
                         selectedPath={selectedFile?.path}
                      />
                   ))}
                </div>
             </div>

             {/* Content Area */}
             <div className="flex-1 flex flex-col min-w-0 bg-white">
                {selectedFile ? (
                   <>
                      {/* File Tab Header */}
                      <div className="flex items-center gap-2 px-4 py-2 border-b border-stone-100 bg-white sticky top-0">
                         {getFileIcon(selectedFile.name)}
                         <span className="text-sm font-medium text-stone-700">{selectedFile.name}</span>
                         <span className="text-xs text-stone-400 ml-auto font-mono">{selectedFile.path}</span>
                      </div>
                      
                      {/* Code/Preview */}
                      <div className="flex-1 overflow-auto custom-scrollbar relative p-0 bg-white">
                         {contentLoading ? (
                            <div className="absolute inset-0 flex items-center justify-center bg-white/80 z-10">
                               <Loader2 size={32} className="animate-spin text-orange-400" />
                            </div>
                         ) : isImage ? (
                            <div className="flex items-center justify-center h-full p-10 bg-stone-50">
                               <img src={fileContent} alt={selectedFile.name} className="max-w-full max-h-full object-contain shadow-lg rounded-lg border border-stone-200" />
                            </div>
                         ) : (
                           selectedFile.name.toLowerCase().endsWith('.md') ? (
                              <div className="prose prose-stone prose-sm max-w-4xl mx-auto p-8">
                                 <ReactMarkdown remarkPlugins={[remarkGfm]}>{fileContent}</ReactMarkdown>
                              </div>
                           ) : (
                             <SyntaxHighlighter 
                                language={selectedFile.language || 'text'} 
                                style={oneLight}
                                customStyle={{ margin: 0, height: '100%', fontSize: '13px', lineHeight: '1.5' }}
                                showLineNumbers={true}
                                lineNumberStyle={{ minWidth: '3em', paddingRight: '1em', color: '#ccc', textAlign: 'right' }}
                             >
                                {fileContent || ''}
                             </SyntaxHighlighter>
                           )
                         )}
                      </div>
                   </>
                ) : (
                   <div className="flex-1 flex flex-col items-center justify-center text-stone-300 gap-4">
                      <div className="p-6 bg-stone-50 rounded-full">
                         <Code size={48} strokeWidth={1} />
                      </div>
                      <p className="font-medium">Select a file to view content</p>
                   </div>
                )}
             </div>
          </div>
       </div>
    </div>
  );
}

function App() {
  const [language, setLanguage] = useState('en');
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedCategory, setSelectedCategory] = useState('All');
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [selectedSkill, setSelectedSkill] = useState(null);

  // Helper to get localized text
  const t = (content) => {
    if (!content) return '';
    if (typeof content === 'string') {
        const simpleDict = {
            'Description': { en: 'Description', zh: '简介' },
            'Use Case': { en: 'Best For', zh: '适用场景' }
        };
        if (simpleDict[content]) return simpleDict[content][language];
        return content;
    }
    return content[language] || content['en'];
  };

  const toggleLanguage = () => {
    setLanguage(prev => {
      const newLang = prev === 'en' ? 'zh' : 'en';
      setSelectedCategory('All');
      return newLang;
    });
  };

  const categories = useMemo(() => {
    return ['All', ...skillsData.map(group => t(group.category))];
  }, [language]);

  const filteredSkills = useMemo(() => {
    let filtered = skillsData.flatMap(group => 
      group.skills.map(skill => ({ 
        ...skill, 
        category: t(group.category),
        description: t(skill.description),
        useCase: t(skill.useCase)
      }))
    );

    if (selectedCategory !== 'All') {
      filtered = filtered.filter(skill => skill.category === selectedCategory);
    }

    if (searchQuery) {
      const query = searchQuery.toLowerCase();
      filtered = filtered.filter(skill => 
        skill.name.toLowerCase().includes(query) || 
        skill.description.toLowerCase().includes(query) ||
        (skill.useCase && skill.useCase.toLowerCase().includes(query))
      );
    }

    return filtered;
  }, [searchQuery, selectedCategory, language]);

  // Handle viewing repo
  // If we have local repo map data, show RepoBrowser. Else fallback to Modal (which tries readme fetching) 
  // - Actually, user wants "whole repo". If we don't have it, we should probably tell them or fallback.
  // We can unify this: "View Source" button -> Opens RepoBrowser.
  
  // NOTE: In this version, repoMap determines if we can show the fancy browser.
  
  return (
    <div className="min-h-screen bg-[#FFFBF7] text-stone-900 font-sans selection:bg-orange-100 selection:text-orange-900">
      
      {/* Modals */}
      {selectedSkill && (
         repoMap[selectedSkill.name] ? (
            <RepoBrowser 
               skill={selectedSkill}
               repoId={repoMap[selectedSkill.name]}
               onClose={() => setSelectedSkill(null)}
            />
         ) : (
            // Fallback for missing repos (or just show simple modal with Warning)
            // For now, let's keep simple modal logic or just a "Not cached" message?
            // Better: Just use the previous logic but inside a modal wrapper, OR
            // Since I replaced the file, I need to restore the Simple Modal as fallback?
            // I'll implement a simple callback here.
            <div className="fixed inset-0 z-[60] flex items-center justify-center p-4 bg-stone-900/50 backdrop-blur-sm" onClick={() => setSelectedSkill(null)}>
               <div className="bg-white p-8 rounded-xl max-w-md w-full text-center" onClick={e => e.stopPropagation()}>
                  <div className="w-16 h-16 bg-orange-100 text-orange-500 rounded-full flex items-center justify-center mx-auto mb-4">
                     <Database size={32} />
                  </div>
                  <h3 className="text-xl font-bold mb-2">Repository Not Cached</h3>
                  <p className="text-stone-600 mb-6">The full repository content hasn't been cached locally yet. You can view it on GitHub.</p>
                  <div className="flex gap-3 justify-center">
                     <button onClick={() => setSelectedSkill(null)} className="px-4 py-2 text-stone-600 font-medium hover:bg-stone-100 rounded-lg">Close</button>
                     <a href={selectedSkill.url} target="_blank" rel="noopener noreferrer" className="px-4 py-2 bg-stone-900 text-white font-medium rounded-lg hover:bg-orange-600 flex items-center gap-2">
                        <ExternalLink size={16} /> Open GitHub
                     </a>
                  </div>
               </div>
            </div>
         )
      )}

      {/* Mobile Sidebar Overlay */}
      {isSidebarOpen && (
        <div 
          className="fixed inset-0 bg-black/20 backdrop-blur-sm z-40 lg:hidden"
          onClick={() => setIsSidebarOpen(false)}
        />
      )}

      {/* Sidebar Navigation */}
      <aside className={`fixed top-0 left-0 z-50 h-full w-72 bg-[#FFFBF7] border-r border-stone-200 transform transition-transform duration-300 ease-in-out lg:translate-x-0 ${isSidebarOpen ? 'translate-x-0' : '-translate-x-full'}`}>
        <div className="h-full flex flex-col">
          <div className="p-6 border-b border-stone-100">
            <div className="flex items-center gap-3">
              <div className="bg-orange-500 text-white p-2.5 rounded-xl shadow-orange-200 shadow-md">
                <Brain size={26} />
              </div>
              <span className="text-xl font-bold text-stone-800">
                {t({ en: 'Edu Skills', zh: '教育技能' })}
              </span>
            </div>
          </div>
          
          <div className="flex-1 overflow-y-auto py-6 px-4 space-y-1 custom-scrollbar">
            <h3 className="px-4 text-xs font-bold text-stone-400 uppercase tracking-widest mb-4">{t({ en: 'Focus Areas', zh: '关注领域' })}</h3>
            {categories.map(category => (
              <button
                key={category}
                onClick={() => {
                  setSelectedCategory(category);
                  setIsSidebarOpen(false);
                }}
                className={`w-full flex items-center gap-3 px-4 py-3.5 text-sm font-medium rounded-xl transition-all duration-200 group
                  ${selectedCategory === category 
                    ? 'bg-orange-50 text-orange-700 shadow-sm ring-1 ring-orange-200' 
                    : 'text-stone-600 hover:bg-white hover:text-stone-900 hover:translate-x-1 hover:shadow-sm'}`}
              >
                {category === 'All' ? <Filter size={18} /> : 
                 CATEGORY_ICONS[category] ? <CategoryIcon icon={CATEGORY_ICONS[category]} size={18} /> : 
                 <Box size={18} />}
                <span className="truncate text-left">{category}</span>
                {category === 'All' && (
                  <span className="ml-auto text-xs bg-stone-100 text-stone-500 py-0.5 px-2 rounded-full group-hover:bg-orange-50 group-hover:text-orange-600 transition-colors">
                    {skillsData.reduce((acc, curr) => acc + curr.skills.length, 0)}
                  </span>
                )}
              </button>
            ))}
          </div>
        </div>
      </aside>

      {/* Main Content */}
      <div className="lg:pl-72 flex flex-col min-h-screen transition-all duration-300">
        
        {/* Header */}
        <header className="sticky top-0 z-30 bg-[#FFFBF7]/90 backdrop-blur-xl border-b border-stone-200 px-4 sm:px-8 py-5">
           <div className="max-w-6xl mx-auto flex items-center justify-between gap-4">
            <button 
              onClick={() => setIsSidebarOpen(true)}
              className="lg:hidden p-2 text-stone-600 hover:bg-white rounded-lg transition-colors"
            >
              <Menu size={24} />
            </button>

            <div className="flex-1 max-w-2xl relative group">
              <div className="absolute inset-y-0 left-0 pl-5 flex items-center pointer-events-none text-stone-400 group-focus-within:text-orange-500 transition-colors">
                <Search size={20} />
              </div>
              <input
                type="text"
                placeholder={t({ en: 'Find educational resources...', zh: '查找教育资源...' })}
                className="block w-full pl-12 pr-4 py-3.5 bg-white border-0 rounded-2xl text-stone-900 shadow-sm ring-1 ring-stone-200 placeholder:text-stone-400 focus:ring-2 focus:ring-orange-500 focus:shadow-lg focus:shadow-orange-100 transition-all duration-200 ease-out"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
              />
            </div>

            <button
                onClick={toggleLanguage}
                className="flex items-center gap-2 px-4 py-2 bg-white border border-stone-200 rounded-2xl text-stone-600 hover:text-orange-600 hover:border-orange-200 hover:shadow-sm transition-all duration-200"
            >
                <Globe size={18} />
                <span className="font-medium text-sm">{language === 'en' ? 'EN' : '中文'}</span>
            </button>
          </div>
        </header>

        {/* Content Area */}
        <main className="flex-1 p-4 sm:p-8 max-w-6xl mx-auto w-full">
           <div className="mb-10">
            <h2 className="text-3xl font-bold text-stone-900 mb-3 tracking-tight">
              {selectedCategory === 'All' ? t({ en: 'All Skills', zh: '所有技能' }) : selectedCategory}
            </h2>
            <p className="text-stone-600 text-base max-w-2xl leading-relaxed">
              {selectedCategory === 'All' 
                ? t({ en: 'Discover a comprehensive curated list of AI resources for education, including MCP servers, specialized LLMs, and Agent frameworks.', zh: '探索全面的教育 AI 资源列表，包含 MCP 服务器、专用大模型和智能体框架。' })
                : t(skillsData.find(c => t(c.category) === selectedCategory)?.description)}
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-6">
            {filteredSkills.map((skill, index) => {
              const Icon = CATEGORY_ICONS[skill.category] || Box;
              const hasRepo = !!repoMap[skill.name];

              return (
                <div key={index} className="group relative bg-white rounded-2xl p-7 shadow-[0_2px_8px_rgba(0,0,0,0.04)] ring-1 ring-stone-100 hover:shadow-[0_20px_40px_-12px_rgba(249,115,22,0.15)] hover:ring-orange-500/30 hover:-translate-y-1 transition-all duration-300 ease-out flex flex-col h-full overflow-hidden cursor-pointer"
                  onClick={() => setSelectedSkill(skill)}
                >
                   {/* Card Top Border Accent */}
                   <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-orange-500 via-amber-500 to-rose-500 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></div>

                  <div className="flex items-start justify-between mb-5">
                    <div className="p-3.5 bg-stone-50 text-stone-600 rounded-2xl group-hover:bg-orange-500 group-hover:text-white transition-all duration-300 group-hover:scale-110 shadow-sm group-hover:shadow-orange-200">
                      <Icon size={26} strokeWidth={1.5} />
                    </div>
                    {hasRepo && (
                        <span className="inline-flex items-center gap-1 px-2 py-1 rounded-md bg-green-50 text-green-700 text-[10px] font-bold uppercase tracking-wider border border-green-100">
                            <Database size={10} /> Local
                        </span>
                    )}
                  </div>

                  <div className="flex-1 mb-6">
                    <h3 className="text-xl font-bold text-stone-900 mb-3 leading-snug group-hover:text-orange-700 transition-colors">
                      {skill.name}
                    </h3>
                    <p className="text-stone-600 text-[15px] leading-relaxed mb-5 line-clamp-3">
                      {skill.description}
                    </p>
                    
                    <div className="flex flex-wrap gap-2">
                       <span className="inline-flex items-center px-3 py-1 rounded-full text-xs font-medium bg-stone-100 text-stone-600 border border-stone-200 group-hover:bg-orange-50 group-hover:border-orange-100 group-hover:text-orange-600 transition-colors">
                          {skill.category.split(' & ')[0]}
                       </span>
                    </div>
                  </div>

                  <div className="pt-6 border-t border-stone-50 mt-auto">
                    {skill.useCase && (
                      <div className="mb-5">
                        <span className="text-xs font-bold text-stone-400 uppercase tracking-wider block mb-1.5">{t({ en: 'Best For', zh: '适用场景' })}</span>
                        <p className="text-sm text-stone-700 font-medium leading-normal bg-stone-50/50 -mx-2 px-2 py-1 rounded-lg">
                          {skill.useCase}
                        </p>
                      </div>
                    )}
                    
                    <button
                        onClick={(e) => {
                            e.stopPropagation();
                            setSelectedSkill(skill);
                        }}
                        className={`flex items-center justify-center w-full py-2.5 px-3 border border-stone-200 text-sm font-semibold rounded-xl transition-all duration-200 gap-2
                             ${hasRepo 
                                ? 'bg-orange-50/50 text-orange-700 border-orange-200 hover:bg-orange-100' 
                                : 'bg-stone-50 text-stone-700 hover:bg-stone-100'}`}
                    >
                        {hasRepo ? <Code size={16} /> : <ExternalLink size={16} />}
                        {hasRepo ? t({ en: 'Browse Code', zh: '浏览代码' }) : t({ en: 'View Details', zh: '查看详情' })}
                    </button>
                  </div>
                </div>
              );
            })}
          </div>
        </main>
      </div>
    </div>
  );
}

const CategoryIcon = ({ icon: Icon, size }) => <Icon size={size} />;

export default App;
