# quiz-master
Simple quiz interface for GK question for students . 

## Main idea
Idea is to use gen ai and claude code to create a real world application 
* Create set of questions on a topic using claude apis and aws lambda
* Store questions/answer using dynamodb database 
* Serve the quiz over a UI 

## Gen ai
Initially i started with prompting and using claude web interface to generate the code.
Now experimenting with claude code to do small feature changes ( e.g add a spinner while api response is coming)
Application is purposely implemented as a plain vanilla javascript application so that complexity is minimal .
 