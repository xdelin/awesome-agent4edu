// Move all the quiz logic here (keeping it the same)
let quizData = [];
let currentQuestion = 0;
let score = 0;
let selectedOption = -1;
let answered = false;
let userAnswers = [];




function loadQuizInternal() {
    console.log('inside loadQuizInternal  ');
    try {
        const jsonText = document.getElementById('jsonData').value;
        if (jsonText) {
            quizData = JSON.parse(jsonText);
        }
        else {
            alert('Please provide valid JSON data in the textarea');
            return;
        }

        console.log('inside loadQuizInternal  '+quizData.length);

        if (!Array.isArray(quizData) || quizData.length === 0) {
            alert('Please provide a valid array of questions');
            return;
        }

        for (let q of quizData) {
            if (!q.question || !q.options || !Array.isArray(q.options) || typeof q.correct !== 'number') {
                alert('Invalid question format. Each question needs: question, options (array), and correct (number)');
                return;
            }
        }
        document.getElementById('quizForm').style.display = 'none';
        document.getElementById('jsonInput').style.display = 'none';
        document.getElementById('quizContainer').style.display = 'block';
        document.getElementById('quizQuestions').style.display = 'block';

        currentQuestion = 0;
        score = 0;
        userAnswers = [];
        updateScore();
        showQuestion();
    } catch (e) {
        alert('Invalid JSON format. Please check your data.');
    }
}

function showQuestion() {
    if (currentQuestion >= quizData.length) {
        showFinalResults();
        return;
    }

    const question = quizData[currentQuestion];
    document.getElementById('questionText').textContent = `${currentQuestion + 1}. ${question.question}`;

    const optionsContainer = document.getElementById('optionsContainer');
    optionsContainer.innerHTML = '';

    question.options.forEach((option, index) => {
        const optionElement = document.createElement('div');
        optionElement.className = 'option';
        optionElement.textContent = option;
        optionElement.onclick = () => selectOption(index);
        optionElement.setAttribute('data-index', index);
        optionsContainer.appendChild(optionElement);
    });

    selectedOption = -1;
    answered = false;
    document.getElementById('submitBtn').disabled = true;
    document.getElementById('nextBtn').style.display = 'none';
    document.getElementById('submitBtn').style.display = 'inline-block';

    updateProgress();
}

function selectOption(index) {
    if (answered) return;

    selectedOption = index;

    document.querySelectorAll('.option').forEach(opt => {
        opt.classList.remove('selected');
    });

    document.querySelector(`[data-index="${index}"]`).classList.add('selected');
    document.getElementById('submitBtn').disabled = false;
}

function submitAnswerInternal() {
    if (selectedOption === -1 || answered) return;

    answered = true;
    const question = quizData[currentQuestion];
    const correctIndex = question.correct;

    userAnswers.push({
        questionIndex: currentQuestion,
        question: question.question,
        selectedOption: selectedOption,
        selectedAnswer: question.options[selectedOption],
        correctOption: correctIndex,
        correctAnswer: question.options[correctIndex],
        isCorrect: selectedOption === correctIndex
    });

    document.querySelectorAll('.option').forEach((opt, index) => {
        opt.style.pointerEvents = 'none';
        if (index === correctIndex) {
            opt.classList.add('correct');
        } else if (index === selectedOption && index !== correctIndex) {
            opt.classList.add('incorrect');
        }
    });

    if (selectedOption === correctIndex) {
        score++;
        updateScore();
    }

    document.getElementById('submitBtn').style.display = 'none';
    document.getElementById('nextBtn').style.display = 'inline-block';
}

function nextQuestionInternal() {
    currentQuestion++;
    showQuestion();
}

function updateScore() {
    document.getElementById('score').textContent = score;
    document.getElementById('total').textContent = quizData.length;
}

function updateProgress() {
    const progress = ((currentQuestion + 1) / quizData.length) * 100;
    document.getElementById('progress').style.width = progress + '%';
}

async function showFinalResults() {
    const percentage = Math.round((score / quizData.length) * 100);
    let message = `Quiz Complete! 🎉<br>Final Score: ${score} / ${quizData.length} (${percentage}%)`;

    if (percentage >= 90) message += '<br>Excellent work! 🏆';
    else if (percentage >= 70) message += '<br>Great job! 👏';
    else if (percentage >= 50) message += '<br>Good effort! 👍';
    else message += '<br>Keep practicing! 💪';

    

    document.getElementById('results').innerHTML = message;
    document.getElementById('results').style.display = 'block';
    document.getElementById('nextBtn').style.display = 'none';
    document.getElementById('restartBtn').style.display = 'inline-block';
    document.getElementById('saveBtn').style.display = 'inline-block';
}

async function saveQuizInternal() {
    console.log('Saving quiz ...');
    const quiz_name = sessionStorage.getItem('quizName') || '';
    const prompt = document.getElementById('prompt').value || '';
    let config = await initializeConfig();
    let apiUrl = config.QUIZ_CREATE_URL;

    const quiz = {
        quiz_name : quiz_name || 'Untitled Quiz',
        prompt : prompt,
        quizData: quizData
    };

    const response = await fetch(apiUrl, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'Authorization': `Basic ${btoa(`${sessionStorage.getItem('quizUser')}:${sessionStorage.getItem('quizPassword')}`)}`
        },
        body: JSON.stringify(quiz)
    });

    if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.json();
    console.log('Results saved:', data);
    return data;
}

function restartQuizInternal() {
    currentQuestion = 0;
    score = 0;
    selectedOption = -1;
    answered = false;
    userAnswers = [];

    document.getElementById('results').style.display = 'none';
    document.getElementById('restartBtn').style.display = 'none';
    document.getElementById('saveBtn').style.display = 'none';

    updateScore();
    showQuestion();
}

function saveQuizInternal() {
    if (quizData.length === 0) {
        console.log('No quiz data available to save.');
        return;
    }

    // Proceed with saving the quiz data
    console.log('Saving quiz data...');
}

async function loadJson() {
    try {
        const response = await fetch('/quiz.json');
        const data = await response.json();
        console.log(data);
        return data;
    } catch (error) {
        console.error('Error loading JSON:', error);
    }
}

window.onload = async function () {
    const sampleData = await loadJson();
    if (document.getElementById('jsonData')) {
        document.getElementById('jsonData').value = JSON.stringify(sampleData, null, 2);
    }
};

// AuthManager class definition
class AuthManager {
    constructor() {
        this.isAuthenticated = false;
        this.users = {
            'admin': 'password123',
            'aditi': 'AditiLovesCurd&Coffee',
            'test': 'test123'
        };
    }

    init() {
        // Check if user is already logged in
        const savedAuth = sessionStorage.getItem('quizAuth');
        if (savedAuth === 'true') {
            this.isAuthenticated = true;
            this.showQuizApp();
        } else {
            this.showLoginPage();
        }

        this.setupEventListeners();
    }

    setupEventListeners() {
        const loginForm = document.getElementById('loginForm');
        const logoutBtn = document.getElementById('logoutBtn');

        if (loginForm) {
            loginForm.addEventListener('submit', (e) => {
                e.preventDefault();
                this.handleLogin();
            });
        }

        if (logoutBtn) {
            logoutBtn.addEventListener('click', () => {
                this.handleLogout();
            });
        }
    }

    handleLogin() {
        const username = document.getElementById('username').value;
        const password = document.getElementById('password').value;
        const errorDiv = document.getElementById('loginError');

        console.log('Attempting login with:', username, password);

        if (this.validateCredentials(username, password)) {
            this.isAuthenticated = true;
            sessionStorage.setItem('quizAuth', 'true');
            sessionStorage.setItem('quizUser', username);
            sessionStorage.setItem('quizPassword', password);
            this.showQuizApp();
            this.clearLoginForm();
            if (errorDiv) errorDiv.style.display = 'none';
        } else {
            if (errorDiv) {
                errorDiv.textContent = 'Invalid username or password';
                errorDiv.style.display = 'block';
            }
        }
    }

    handleLogout() {
        this.isAuthenticated = false;
        sessionStorage.removeItem('quizAuth');
        sessionStorage.removeItem('quizUser');
        sessionStorage.removeItem('quizPassword');
        console.log('User logged out');
        this.showLoginPage();
        this.clearLoginForm();
    }

    validateCredentials(username, password) {
        return this.users[username] && this.users[username] === password;
    }

    showLoginPage() {
        const loginPage = document.getElementById('loginPage');
        const quizApp = document.getElementById('quizApp');
        
        if (loginPage) loginPage.style.display = 'block';
        if (quizApp) quizApp.style.display = 'none';
    }

    showQuizApp() {
        const loginPage = document.getElementById('loginPage');
        const quizApp = document.getElementById('quizApp');
        
        if (loginPage) loginPage.style.display = 'none';
        if (quizApp) quizApp.style.display = 'block';
    }

    clearLoginForm() {
        const username = document.getElementById('username');
        const password = document.getElementById('password');
        
        if (username) username.value = '';
        if (password) password.value = '';
    }

    checkAuth() {
        return this.isAuthenticated;
    }
}

// Global variable to hold the auth manager instance
let authManager;

// Enhanced navigation with auth check
function navigateWithAuth(callback) {
    if (authManager && authManager.checkAuth()) {
        callback();
    } else {
        alert('Please log in to access this feature');
        if (authManager) {
            authManager.showLoginPage();
        }
    }
}

// Wait for DOM to be ready - THIS IS WHERE WE INITIALIZE EVERYTHING
document.addEventListener('DOMContentLoaded', function() {
    // Initialize authentication manager AFTER DOM is loaded
    authManager = new AuthManager();
    authManager.init();

    // Menu toggle functionality
    const menuToggle = document.getElementById('menuToggle');
    const sideMenu = document.getElementById('sideMenu');
    const overlay = document.getElementById('overlay');

    if (menuToggle && sideMenu && overlay) {
        menuToggle.addEventListener('click', function() {
            navigateWithAuth(() => {
                console.log('Menu toggle clicked');
                sideMenu.classList.toggle('active');
                overlay.style.display = sideMenu.classList.contains('active') ? 'block' : 'none';
            });
        });

        overlay.addEventListener('click', function() {
            sideMenu.classList.remove('active');
            overlay.style.display = 'none';
        });
    }

    // Navigation links with auth protection
    const createQuizLink = document.getElementById('createQuizLink');
    const homeLink = document.getElementById('homeLink');
    const quizForm = document.getElementById('quizForm');
    const quizQuestions = document.getElementById('quizQuestions');

    if (createQuizLink) {
        createQuizLink.addEventListener('click', function(e) {
            e.preventDefault();
            navigateWithAuth(() => {
                if (quizForm) quizForm.style.display = 'block';
                if (quizQuestions) quizQuestions.style.display = 'none';
                if (sideMenu) sideMenu.classList.remove('open');
                if (overlay) overlay.style.display = 'none';
            });
        });
    }

    if (homeLink) {
        homeLink.addEventListener('click', function(e) {
            e.preventDefault();
            navigateWithAuth(() => {
                if (quizForm) quizForm.style.display = 'none';
                if (quizQuestions) quizQuestions.style.display = 'block';
                if (sideMenu) sideMenu.classList.remove('open');
                if (overlay) overlay.style.display = 'none';
            });
        });
    }

    // Handle quiz form submission
    if (quizForm) {
        quizForm.addEventListener('submit', async (e) => {
            e.preventDefault();
            const promptInput = document.getElementById('prompt');
            if (promptInput) {
                // call api with auth and prompt
                config = await initializeConfig();
                const prompt = promptInput.value;
                if (prompt.trim()) {
                    console.log('Quiz prompt:', prompt);
                    sessionStorage.setItem('quizPrompt', prompt);
                    
                    // Show spinner before API call
                    showSpinner();
                    
                    try {
                        let username = sessionStorage.getItem('quizUser');
                        let password = sessionStorage.getItem('quizPassword');
                        const apiUrl = config.QUIZ_FETCH_URL;
                        console.log('API URL:', apiUrl);
                        
                        const quizData = await callBackendAPI(apiUrl, username, password, prompt);
                        console.log('Quiz data fetched:', quizData);
                        
                        if (quizData && quizData["data"]) {
                            if (document.getElementById('jsonData')) {
                                console.log('quizData set to html field - jsonData');
                                document.getElementById('jsonData').value = JSON.stringify(quizData["data"], null, 2);
                            }
                            
                            loadQuizInternal();
                        } else {
                            alert('Failed to fetch quiz data. Please try again.');
                        }
                    } catch (error) {
                        console.error('Error generating quiz:', error);
                        alert('Error generating quiz. Please try again.');
                    } finally {
                        // Always hide spinner, whether success or error
                        hideSpinner();
                    }
                    
                } else {
                    alert('Please enter a prompt for your quiz.');
                }
            }
        });
    }
});

// Override existing functions to include auth checks
const originalFunctions = {
    loadQuiz: window.loadQuiz,
    submitAnswer: window.submitAnswer,
    nextQuestion: window.nextQuestion,
    restartQuiz: window.restartQuiz,
    saveQuiz : window.saveQuiz
};

window.loadQuiz = function() {
    if (authManager && authManager.checkAuth()) {
        loadQuizInternal();
    } else {
        alert('Please log in to load a quiz');
    }
};

window.submitAnswer = function() {
    if (authManager && authManager.checkAuth()) {
        submitAnswerInternal();
    }
};

window.nextQuestion = function() {
    if (authManager && authManager.checkAuth()) {
        nextQuestionInternal();
    }
};

window.restartQuiz = function() {
    if (authManager && authManager.checkAuth()) {
        restartQuizInternal();
    }
};

window.restartQuiz = function() {
    if (authManager && authManager.checkAuth()) {
        restartQuizInternal();
    }
};

window.saveQuiz = function() {
    if (authManager && authManager.checkAuth()) {
        saveQuizInternal();
    }
};

// Spinner control functions
function showSpinner() {
    const spinnerOverlay = document.getElementById('spinnerOverlay');
    if (spinnerOverlay) {
        spinnerOverlay.style.display = 'flex';
    }
}

function hideSpinner() {
    const spinnerOverlay = document.getElementById('spinnerOverlay');
    if (spinnerOverlay) {
        spinnerOverlay.style.display = 'none';
    }
}

async function callBackendAPI(url, username, password, promptText) {
  // Create basic auth header
  const credentials = btoa(`${username}:${password}`);
  /*   'Authorization': `Basic ${credentials}` */
  try {
    const response = await fetch(url, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Basic ${credentials}`
      },
      body: JSON.stringify({
        prompt: promptText
      })
    });
    
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    
    const data = await response.json();
    return data;
    
  } catch (error) {
    console.error('API call failed:', error);
    throw error;
  }
}


// config.js
let appConfig = null;
let isInitialized = false;

async function initializeConfig() {
  if (isInitialized) return appConfig;
  
  try {
    const response = await fetch('./config.json');
    appConfig = await response.json();
  } catch (error) {
    appConfig = {
      quizCreateUrl: 'http://localhost:3000',
      quizFetchUrl: 'http://localhost:3000',
      environment: 'development'
    };
  }
  isInitialized = true;
  return appConfig;
}

