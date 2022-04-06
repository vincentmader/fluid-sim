
run:
	rm ./data/*/* > /dev/null
	cargo run
	rm ./plots/*/* > /dev/null
	./src/visualize/main.py

calc:
	rm ./data/*/* > /dev/null
	cargo run

vis:
	rm ./plots/*/* > /dev/null
	./src/visualize/main.py

