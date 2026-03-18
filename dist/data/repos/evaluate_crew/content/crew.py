from crewai import Agent, Crew, Process, Task
from crewai.project import CrewBase, agent, crew, task

# If you want to run a snippet of code before or after the crew starts, 
# you can use the @before_kickoff and @after_kickoff decorators
# https://docs.crewai.com/concepts/crews#example-crew-class-with-decorators
from typing import List, Optional
from pydantic import BaseModel, Field
from datetime import datetime, timedelta


class quiz(BaseModel):
    score : int = Field(..., description="Number of correct answer out of 10 questions", ge=0,e=10,)
    performance : str = Field(..., description="Don't provide result just provide Overall performance, including areas of strength and areas needing improvement")
    
class AnalysisOutput(BaseModel):
    output : List[quiz] = Field(..., description= "Perfromace score and over all topics to learn")



@CrewBase
class Evaluate():
	"""Evaluate crew"""

	# Learn more about YAML configuration files here:
	# Agents: https://docs.crewai.com/concepts/agents#yaml-configuration-recommended
	# Tasks: https://docs.crewai.com/concepts/tasks#yaml-configuration-recommended
	agents_config = 'config/agents.yaml'
	tasks_config = 'config/tasks.yaml'

	# If you would like to add tools to your agents, you can learn more about it here:
	# https://docs.crewai.com/concepts/agents#agent-tools
	@agent
	def Quiz_bot(self) -> Agent:
		return Agent(
			config=self.agents_config['Quiz_bot'],
			allow_delegation = False,
			verbose=True
		)

	@agent
	def analyze_bot(self) -> Agent:
		return Agent(
			config=self.agents_config['analyze_bot'],
			verbose=True,
		)

	@agent
	def quality_agent(self) -> Agent:
		return Agent(
			config=self.agents_config['quality_agent'],
   			allow_delegation = True,
			verbose=True,
		)
  
  
	

	# To learn more about structured task outputs, 
	# task dependencies, and task callbacks, check out the documentation:
	# https://docs.crewai.com/concepts/tasks#overview-of-a-task
	@task
	def quiz_task(self) -> Task:
		return Task(
			config=self.tasks_config['quiz_task'],
			memory = True,
			human_input = True,
			output_pydantic = AnalysisOutput,
   
			
		)

	@task
	def analyze_task(self) -> Task:
		return Task(
			config=self.tasks_config['analyze_task'],
			
			
			
		)
	
	@task
	def quality_task(self) -> Task:
		return Task(
			config=self.tasks_config['quality_task'],
			output_file = "report.md",
			
		)
  
	

	@crew
	def crew(self) -> Crew:
		"""Creates the Evaluate crew"""
		# To learn how to add knowledge sources to your crew, check out the documentation:
		# https://docs.crewai.com/concepts/knowledge#what-is-knowledge

		return Crew(
			agents=self.agents, # Automatically created by the @agent decorator
			tasks=self.tasks, # Automatically created by the @task decorator
			process=Process.sequential,
			verbose=True,
			# process=Process.hierarchical, # In case you wanna use that instead https://docs.crewai.com/how-to/Hierarchical/
		)
