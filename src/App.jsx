
import { useState, useEffect, useMemo } from 'react';
import skillsData from './data/skills.json';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { Search, ExternalLink, BookOpen, Code, PenTool, Database, Layout, Brain, Calculator, Briefcase, Box, Menu, X, Filter, Globe, Loader2, Info } from 'lucide-react';

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

function SkillDetailModal({ skill, onClose, t }) {
  const [readme, setReadme] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    if (!skill) return;

    const fetchRepoDetails = async () => {
      setLoading(true);
      setError(null);
      setReadme(null);

      try {
        // Simple extraction of owner/repo from URL
        // Expects https://github.com/owner/repo
        const match = skill.url.match(/github\.com\/([^/]+)\/([^/]+)/);
        if (!match) {
          // If not github, maybe we can't fetch readme easily
          throw new Error("Preview available only for GitHub repositories.");
        }
        
        const [, owner, repo] = match;
        
        // Fetch README via GitHub API
        const response = await fetch(`https://api.github.com/repos/${owner}/${repo}/readme`);
        
        if (response.status === 403) {
           throw new Error("GitHub API Rate Limit Exceeded. Please view on GitHub directly.");
        }
        
        if (!response.ok) {
           throw new Error("Could not load README.");
        }

        const data = await response.json();
        // Decode Base64 securely for UTF-8
        try {
            const content = decodeURIComponent(escape(window.atob(data.content.replace(/\s/g, ''))));
            setReadme(content);
        } catch (e) {
            setReadme(window.atob(data.content.replace(/\s/g, '')));
        }
        
      } catch (err) {
        setError(err.message || "Failed to load repository details");
      } finally {
        setLoading(false);
      }
    };

    fetchRepoDetails();
  }, [skill]);

  if (!skill) return null;

  return (
    <div className="fixed inset-0 z-[60] flex items-center justify-center p-4 sm:p-6" role="dialog" aria-modal="true">
      {/* Backdrop */}
      <div 
        className="fixed inset-0 bg-stone-900/50 backdrop-blur-sm transition-opacity" 
        onClick={onClose}
      ></div>

      {/* Modal Content */}
      <div className="relative bg-white rounded-2xl shadow-2xl w-full max-w-4xl max-h-[90vh] flex flex-col overflow-hidden animate-in fade-in zoom-in-95 duration-200">
        
        {/* Header */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-stone-100 bg-stone-50/50">
          <div className="flex items-center gap-3">
             <div className="bg-orange-100/50 p-2 rounded-lg text-orange-600">
                {CATEGORY_ICONS[skill.category] ? <CategoryIcon icon={CATEGORY_ICONS[skill.category]} size={20} /> : <Box size={20} />}
             </div>
             <div>
                <h3 className="text-lg font-bold text-stone-900 leading-snug">{skill.name}</h3>
                <p className="text-xs text-stone-500 font-medium">{skill.category}</p>
             </div>
          </div>
          <button 
            onClick={onClose}
            className="p-2 text-stone-400 hover:text-stone-600 hover:bg-stone-100 rounded-full transition-colors"
          >
            <X size={20} />
          </button>
        </div>

        {/* Scrollable Content */}
        <div className="flex-1 overflow-y-auto p-6 custom-scrollbar bg-white">
           
           {/* Basic Info Block */}
           <div className="mb-8 p-5 bg-orange-50/30 border border-orange-100/50 rounded-xl">
              <h4 className="text-sm font-bold text-stone-800 uppercase tracking-wide mb-2 flex items-center gap-2">
                 <Info size={14} className="text-orange-500"/> 
                 {t('Description')}
              </h4>
              <p className="text-stone-700 text-sm leading-relaxed mb-4">{skill.description}</p>
              
              {skill.useCase && (
                  <>
                     <h4 className="text-sm font-bold text-stone-800 uppercase tracking-wide mb-2 mt-4">{t('Use Case')}</h4>
                     <p className="text-stone-700 text-sm leading-relaxed">{skill.useCase}</p>
                  </>
              )}
           </div>

           {/* README Loader/Content */}
           <div className="prose prose-stone prose-sm max-w-none prose-headings:font-bold prose-h1:text-2xl prose-h2:text-xl prose-a:text-orange-600 hover:prose-a:text-orange-700 prose-img:rounded-xl">
              {loading && (
                 <div className="flex flex-col items-center justify-center py-12 text-stone-400">
                    <Loader2 size={32} className="animate-spin mb-3 text-orange-400" />
                    <p className="text-sm">Loading repository details...</p>
                 </div>
              )}

              {error && (
                 <div className="text-center py-10 bg-stone-50 rounded-xl border border-stone-100 border-dashed">
                    <p className="text-stone-500 mb-4">{error}</p>
                    <a href={skill.url} target="_blank" rel="noopener noreferrer" className="inline-flex items-center gap-2 px-4 py-2 bg-white border border-stone-200 shadow-sm rounded-lg text-sm font-medium hover:bg-stone-50 transition-colors">
                       <ExternalLink size={14} />
                       View on GitHub
                    </a>
                 </div>
              )}

              {!loading && !error && readme && (
                  <ReactMarkdown remarkPlugins={[remarkGfm]}>
                      {readme}
                  </ReactMarkdown>
              )}
           </div>
        </div>

        {/* Footer */}
        <div className="p-4 bg-stone-50 border-t border-stone-100 flex justify-end gap-3 z-10">
           <button 
             onClick={onClose}
             className="px-4 py-2 text-sm font-medium text-stone-600 hover:bg-stone-200/50 rounded-lg transition-colors"
           >
             Close
           </button>
           <a 
             href={skill.url} 
             target="_blank" 
             rel="noopener noreferrer"
             className="px-4 py-2 text-sm font-semibold text-white bg-stone-900 hover:bg-orange-600 rounded-lg shadow-sm hover:shadow-md transition-all flex items-center gap-2"
           >
             <Code size={16} />
             View Source
           </a>
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
      setSelectedCategory('All'); // Reset category to avoid mismatch
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

  // Labels for static UI elements
  const uiLabels = {
    allSkills: { en: 'All Skills', zh: '所有技能' },
    allSkillsDesc: { en: 'Discover a curated collection of AI skills for educators, students, and researchers.', zh: '探索为教育工作者、学生和研究人员精选的 AI 技能集合。' },
    focusAreas: { en: 'Focus Areas', zh: '关注领域' },
    learningHub: { en: 'Learning Hub', zh: '学习中心' },
    learningHubDesc: { en: 'Explore AI-powered tools to enhance teaching and learning.', zh: '探索增强教学和学习体验的 AI 工具。' },
    findResources: { en: 'Find educational resources...', zh: '查找教育资源...' },
    bestFor: { en: 'Best For', zh: '适用场景' },
    sourceCode: { en: 'Source Code', zh: '源代码' },
    viewDetails: { en: 'View Details', zh: '查看详情' },
    noSkillsFound: { en: 'No skills found', zh: '未找到相关技能' },
    tryAdjusting: { en: 'Try adjusting your search terms or category.', zh: '尝试调整搜索关键词或类别。' },
    appName: { en: 'Edu Skills', zh: '教育技能' }
  };

  return (
    <div className="min-h-screen bg-[#FFFBF7] text-stone-900 font-sans selection:bg-orange-100 selection:text-orange-900">
      
      {/* Skill Detail Modal */}
      {selectedSkill && (
        <SkillDetailModal 
            skill={selectedSkill} 
            onClose={() => setSelectedSkill(null)} 
            t={t}
        />
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
                {t(uiLabels.appName)}
              </span>
            </div>
          </div>
          
          <div className="flex-1 overflow-y-auto py-6 px-4 space-y-1 custom-scrollbar">
            <h3 className="px-4 text-xs font-bold text-stone-400 uppercase tracking-widest mb-4">{t(uiLabels.focusAreas)}</h3>
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

          <div className="p-4 border-t border-stone-100">
            <div className="bg-gradient-to-br from-orange-50 to-amber-50 rounded-xl p-5 border border-orange-100">
               <h4 className="text-sm font-semibold text-orange-900 mb-1">{t(uiLabels.learningHub)}</h4>
               <p className="text-xs text-orange-800/80 mb-0">
                 {t(uiLabels.learningHubDesc)}
               </p>
            </div>
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
                placeholder={t(uiLabels.findResources)}
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
              {selectedCategory === 'All' ? t(uiLabels.allSkills) : selectedCategory}
            </h2>
            <p className="text-stone-600 text-base max-w-2xl leading-relaxed">
              {selectedCategory === 'All' 
                ? t(uiLabels.allSkillsDesc)
                : t(skillsData.find(c => t(c.category) === selectedCategory)?.description)}
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-6">
            {filteredSkills.map((skill, index) => {
              const Icon = CATEGORY_ICONS[skill.category] || Box;
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
                        <span className="text-xs font-bold text-stone-400 uppercase tracking-wider block mb-1.5">{t(uiLabels.bestFor)}</span>
                        <p className="text-sm text-stone-700 font-medium leading-normal bg-stone-50/50 -mx-2 px-2 py-1 rounded-lg">
                          {skill.useCase}
                        </p>
                      </div>
                    )}
                    
                    <div className="grid grid-cols-2 gap-3">
                         <button
                            onClick={(e) => {
                                e.stopPropagation();
                                setSelectedSkill(skill);
                            }}
                            className="flex items-center justify-center w-full py-2.5 px-3 bg-stone-50 border border-stone-200 text-stone-700 text-sm font-semibold rounded-xl hover:bg-stone-100 transition-all duration-200 gap-2"
                        >
                            <BookOpen size={16} />
                            {t(uiLabels.viewDetails)}
                        </button>

                        <a
                        href={skill.url}
                        target="_blank" 
                        rel="noopener noreferrer"
                        onClick={(e) => e.stopPropagation()}
                        className="flex items-center justify-center w-full py-2.5 px-3 bg-white border-2 border-stone-100 text-stone-700 text-sm font-semibold rounded-xl hover:bg-stone-800 hover:border-stone-800 hover:text-white transition-all duration-200 gap-2"
                        >
                        <Code size={16} />
                        {t(uiLabels.sourceCode)}
                        </a>
                    </div>
                  </div>
                </div>
              );
            })}
          </div>

          {filteredSkills.length === 0 && (
            <div className="text-center py-20">
              <div className="inline-flex items-center justify-center w-16 h-16 rounded-full bg-slate-50 mb-4">
                <Search size={32} className="text-slate-300" />
              </div>
              <h3 className="text-lg font-medium text-slate-900">{t(uiLabels.noSkillsFound)}</h3>
              <p className="text-slate-500">{t(uiLabels.tryAdjusting)}</p>
            </div>
          )}
        </main>
      </div>
    </div>
  );
}

// Helper component to render icons properly
const CategoryIcon = ({ icon: Icon, size }) => <Icon size={size} />;

export default App;
